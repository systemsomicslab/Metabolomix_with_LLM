# 分子種ごとの図・分子種 PCA・ローディング図・図の保存の一本化 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** ARF 経路に分子種ごとの図（割合 % / PeakHeight）と分子種を選んで回す PCA を足し、3 経路の PCA のローディング図と、図の保存ツールの一本化（`save_figure`）を実装する。

**Architecture:** 項目（クラス・分子種名）の解決を `arf/selection.py` に切り出して群別強度・分子種の図・分子種 PCA で共有する。描画は MCP 非依存の `plots/{species,pca_scores,pca_loadings}.py`、MCP ツールは `arf/species_tools.py`・`tools/pca_loadings_tools.py`・`tools/reports.py`。既存の 4 つの保存ツールは `save_figure` に統合したことを確かめてから消す。

**Tech Stack:** Python 3.14（`C:/Python314/python.exe`）、FastMCP、numpy、scikit-learn、matplotlib、pytest。

**Spec:** `docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md`

## Global Constraints

- 描画 payload の形式名: `lipidmix.species_intensity.v1`、`lipidmix.pca_loadings.v1`（外部契約の旧名規則）。新しい環境変数は作らない。
- 全ツールに `structured_output=False`。JSON は `metabolomix.core.serialization.json_payload()` で返す。float は丸める。
- 描画ツールの既定は画像（`[説明文, Image(png)]`）、`output="payload"` か `LIPIDMIX_PLOT_OUTPUT=payload` で payload（`plots/render.py resolve_plot_output()`）。画像は `figure_to_png()`（dpi 100）。
- 前提の状態がないときは `mcp_errors.missing_state(...)`、入力の誤りは `{"status": "error", "message": ...}`。失敗した描画呼び出しの後は、その種類の直近結果を `None` にする。
- 分子種の図のパネルは 1 回 40 枚まで。ローディングの全件表示は 60 特徴量まで。`pcs` は 1〜3 個。
- `save_figure` の `kind`: `pca` / `volcano` / `eic` / `group_intensity` / `species` / `pca_loadings`。`source` / `result_id` は `pca` / `volcano` だけ。`annotations` は `readOnlyHint=False, destructiveHint=False, idempotentHint=True`。
- `analysis/pca.py` の既存 `run_pca` は変えない（`arf/reader.py` の再エクスポートとテストのモック）。
- 図の文字は ASCII（matplotlib の既定フォントは日本語の字形を持たない）。説明文（戻り値の文字列）は日本語でよい。
- コミット時に pre-commit が全テストを走らせる（約 3〜4 分）。`git commit` はバックグラウンドで実行し、出力をファイルへ落として `FAILED` を確認する。実行中は作業ツリーを編集しない。
- コミットメッセージの末尾に `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`。
- 作業ブランチは `feat/species-plot-and-pca-loadings`（作成済み）。
- 手順中の `<scratchpad>` は、そのセッションの scratchpad ディレクトリ（リポジトリの外）を指す。コミットのログはそこへ落とす。

## ツール数が変わるタスクで必ず直す 7 か所（以下「台帳 7 か所」）

ツールを足す・消すタスク（Task 2 / 4 / 7 / 10）は、同じコミットの中で次をすべて直す。どれかを忘れると pre-commit の腐敗防止テストで落ちる。

1. `tests/test_server_registration.py` の `EXPECTED_TOOLS`（と経緯コメント 1 行）
2. `tests/test_tool_annotations.py` の対応表（新ツールの `annotations`）
3. `tests/test_workflow_docs.py` の `IN_SCOPE`（文書ごとの対象ツール）と末尾の合計リテラル（現在 `53`）
4. `docs/workflow/<文書>.md` の `## <ツール名>` 節（前提・状態変更・呼び出し連鎖。連鎖の関数は実在すること）
5. `docs/workflow/index.md` の目次の件数・「合計 N ツール」・「登録ツール総数は N（N + 21）」
6. `USAGE.md` の表の行と見出しの「全Nツール」
7. `CLAUDE.md` 冒頭の「ツール N・リソース 5・リソーステンプレート 3」

ツール数の推移: 74 →（Task 2）71 →（Task 4）72 →（Task 7）73 →（Task 10）74。対象範囲（`IN_SCOPE` の合計）の推移: 53 → 50 → 51 → 52 → 53。

## Review Focus

1. **同じ名前で付加イオンが違うスポット**（例: DG 34:1 の [M+NH4]+ と [M+Na]+）は別パネルになり、見出しの付加イオンと ID で区別できること（Task 3 / Task 4 のテスト）。
2. **ある群の試料がすべて低信頼**のとき、その群の棒・平均線が描かれないだけで落ちないこと（Task 3 のテスト）。
3. **`share_basis` が描く分子種を含まない**（例: PG を描き分母は DGDG だけ）とき、計算は続けるが説明文と payload の caveats に「分母に含まれない分子種がある」と出ること（Task 4 のテスト）。
4. **`normalize="total"` で選んだ分子種の合計が 0 の試料**があるとき、その試料名を挙げたエラーになること（Task 7 のテスト）。
5. **ローディングを持たない古い ARF の PCA 結果**（この変更前に作られたもの）や、**存在しない主成分**（`pcs=[4]` など）を指定したとき、落ちずに分かるエラーになること（Task 10 のテスト）。

---

### Task 1: 項目解決の切り出し（`arf/selection.py`）と群別強度の置き換え

**Files:**
- Create: `metabolomix/arf/selection.py`
- Modify: `metabolomix/arf/tools.py`（`arf_plot_group_intensity` の本体）
- Modify: `docs/workflow/arf.md`（`## arf_plot_group_intensity` の連鎖）
- Test: `tests/test_arf_selection.py`（新規）、既存 `tests/test_group_intensity_tools.py`

**Interfaces:**
- Produces:
  - `selection.Selection`（dataclass）: `items: list[dict]`（`group_intensity.resolve_items` の戻り値の 1 つ目）、`excluded: dict[str, list[str]]`、`groups: list[dict]`（`{"label", "samples"}`）、`low_reliability: frozenset[str]`、`rows_by_spot: dict[int, list[dict]]`（`arf_exclude` した試料を除いた行）、`msms: dict[int, bool]`、`catalog_by_id: dict[int, dict]`、`caveats: list[str]`
  - `selection.build_selection(arf_state, arf2_path, *, items, groups, low_reliability_samples=None, apply_curation=True, exclude_auto_likely_wrong=False, standard_samples=None) -> Selection`。入力の誤りは `ValueError`、判断の記録ファイルが読めなければ `curation.flags.FlagFileError` をそのまま投げる。
  - `selection.expand_spots(selection, *, require_msms=False) -> tuple[list[dict], list[str]]`。1 要素 = `{"spot_id", "name", "label", "ontology", "item", "adduct", "mz", "rt", "msms"}`。2 つ目は MS/MS が無くて外したスポットのラベル `"#<id> <name>"`（昇順）。
  - `selection.display_name(name: str) -> str`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_arf_selection.py`:

```python
"""arf/selection.py: 項目の解決とスポット単位への展開。fixture はテストが作る。"""
from __future__ import annotations

import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state
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
        _group([10, 5000, 6000, 2000, 2500, 10]),            # 0 PG 34:1（MS/MS あり）
        _group([10, 3000, 3500, 4000, 4200, 10]),            # 1 PG 35:1（MS/MS なし）
        _group([0, 900, 950, 800, 700, 0]),                  # 2 DG 34:1 [M+NH4]+
        _group([0, 300, 320, 280, 260, 0]),                  # 3 DG 34:1 [M+Na]+
    ])
    arf2 = write_arf2(tmp_path / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=0, name="PG 34:1|PG 16:0_18:1", ontology="PG", adduct="[M-H]-",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=1, name="PG 35:1|PG 16:0_19:1", ontology="PG", adduct="[M-H]-"),
        arf2_spot_raw(spot_id=2, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+NH4]+",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=3, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+Na]+",
                      matches=[match_result()]),
    ])
    session_state.session = session_state.AnalysisSession()
    with open(arf, "rb") as handle:
        session_state.session.arf.features = arf_reader.deserialize(handle)
    session_state.session.arf.current_file_path = str(arf)
    session_state.session.arf.class_index = None
    return {"arf": arf, "arf2": arf2}


def _select(loaded, **kw):
    from metabolomix.arf import selection
    kw.setdefault("groups", ["ctrl", "ko"])
    return selection.build_selection(session_state.session.arf, loaded["arf2"], **kw)


def test_build_selection_resolves_items_groups_and_msms(loaded):
    sel = _select(loaded, items=["PG", "DG"])
    assert [s["spot_id"] for s in sel.items[0]["spots"]] == [0, 1]
    assert [g["label"] for g in sel.groups] == ["ctrl", "ko"]
    assert sel.msms[0] is True and sel.msms.get(1, False) is False
    assert sel.catalog_by_id[2]["AdductType"] == "[M+NH4]+"


def test_expand_spots_orders_by_item_then_name_and_keeps_adduct_duplicates(loaded):
    from metabolomix.arf import selection
    spots, dropped = selection.expand_spots(_select(loaded, items=["DG", "PG"]))
    assert [s["spot_id"] for s in spots] == [2, 3, 0, 1]
    assert [s["adduct"] for s in spots[:2]] == ["[M+NH4]+", "[M+Na]+"]
    assert spots[0]["label"] == "DG 16:0_18:1" and spots[0]["item"] == "DG"
    assert dropped == []


def test_expand_spots_require_msms_drops_unsupported(loaded):
    from metabolomix.arf import selection
    spots, dropped = selection.expand_spots(_select(loaded, items=["PG"]), require_msms=True)
    assert [s["spot_id"] for s in spots] == [0]
    assert dropped == ["#1 PG 35:1|PG 16:0_19:1"]


def test_expand_spots_counts_a_spot_once_across_items(loaded):
    from metabolomix.arf import selection
    spots, _ = selection.expand_spots(_select(loaded, items=["PG", "PG 34:1"]))
    assert [s["spot_id"] for s in spots] == [0, 1]


def test_display_name_takes_the_last_candidate():
    from metabolomix.arf.selection import display_name
    assert display_name("PG 34:1|PG 16:0_18:1") == "PG 16:0_18:1"
    assert display_name("no MS2: PC 34:1") == "no MS2: PC 34:1"
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_arf_selection.py -q`
Expected: FAIL（`ModuleNotFoundError: No module named 'metabolomix.arf.selection'`）

- [ ] **Step 3: `metabolomix/arf/selection.py` を書く**

```python
"""ARF の項目指定（クラス・分子種名）を、描画・PCA が使うスポットの集合に解決する。

群別強度（arf_plot_group_intensity）・分子種ごとの図（arf_plot_species）・分子種 PCA
（arf_pca_species）が同じ規則で分子種を選ぶための共通部品。項目の文字列の解決そのものは
`plots/group_intensity.py resolve_items()`。ここは ARF の行・`.arf2`・キュレーションの判断と
最新レビュー・標準液の判定を集めて渡す層。MCP には依存しない（session の ArfState は引数で受ける）。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §3。
"""
from __future__ import annotations

from dataclasses import dataclass, field

from metabolomix.arf import reader as arf_reader
from metabolomix.arf2 import reader as arf2_reader
from metabolomix.arf2.match_results import load_spot_annotations, name_prefix
from metabolomix.curation import apply as curation_apply
from metabolomix.curation import flags as curation_flags
from metabolomix.curation import suggest as curation_suggest
from metabolomix.msdial import sample_factors
from metabolomix.plots import group_intensity as gi


@dataclass
class Selection:
    items: list
    excluded: dict
    groups: list
    low_reliability: frozenset
    rows_by_spot: dict
    msms: dict
    catalog_by_id: dict
    caveats: list = field(default_factory=list)


def display_name(name: str) -> str:
    """`PG 34:1|PG 16:0_18:1` → `PG 16:0_18:1`（鎖まで決まった最後の候補）。接頭辞は残す。"""
    return (name or "").split("|")[-1].strip()


def build_selection(arf_state, arf2_path, *, items, groups, low_reliability_samples=None,
                    apply_curation=True, exclude_auto_likely_wrong=False, standard_samples=None):
    """項目・群・除外を解決する。入力の誤りは ValueError、記録ファイルの破損は FlagFileError。"""
    excluded_samples = set(arf_state.excluded_samples or ())
    rows_by_spot = {}
    all_rows_by_spot = {}    # 標準液の判定用: arf_exclude 済みの試料の行も含める
    for feature in arf_state.features:
        rows = [r for r in (arf_reader.alignment_feature_row(raw) for raw in feature["AlignedPeakProperties"]) if r]
        all_rows_by_spot[feature["MasterAlignmentID"]] = rows
        rows_by_spot[feature["MasterAlignmentID"]] = [
            r for r in rows if r.get("file_name") not in excluded_samples]
    names = [n for n in sample_factors.arf_sample_names(arf_state.features) if n not in excluded_samples]
    facets = sample_factors.build_sample_facets(names, arf_state.class_index)
    resolved_groups, caveats = gi.resolve_groups(groups, facets)
    plotted = {n for g in resolved_groups for n in g["samples"]}
    low = gi.resolve_sample_specs(low_reliability_samples or [], facets)

    catalog = arf2_reader.load_catalog(str(arf2_path))
    curation = {}
    if apply_curation:
        flag_set = curation_apply.flags_for_arf2(arf2_path)
        curation = {
            "wrong": {int(s) for s in flag_set["wrong"]},
            "redundant": {int(s) for s in flag_set["redundant"]},
            "assign": {int(s): {"name": row.get("name") or "", "ontology": row.get("ontology") or ""}
                       for s, row in flag_set["assign"].items()},
        }
        if flag_set["orphaned"]:
            caveats.append(curation_flags.orphaned_warning(flag_set["orphaned"]))
    if exclude_auto_likely_wrong:
        review = curation_suggest.latest_review(
            arf2_path, curation_flags.alignment_key(arf2_path)["alignment_sha256"])
        if review is None:
            caveats.append("このアラインメントの curation_review のレビューが無いため、自動判定 likely_wrong は除いていません。")
        else:
            curation["auto_likely_wrong"] = {
                int(s["spot_id"]) for s in review["spots"] if s.get("verdict") == "likely_wrong"}

    standard = set()
    if standard_samples:
        standard = gi.resolve_sample_specs(standard_samples, sample_factors.build_sample_facets(
            sample_factors.arf_sample_names(arf_state.features), arf_state.class_index))
    standard_only = gi.standard_only_spots(all_rows_by_spot, standard, plotted) if standard else frozenset()
    manual = {int(s) for s in (arf_state.excluded_spots or ())}    # arf_exclude で除いたスポット
    resolved_items, excluded = gi.resolve_items(
        items, catalog, curation=curation, standard_only=standard_only, manual=manual)

    # MS/MS の裏付け: 照合結果に MS/MS があり、かつ .arf2 の Name が `no MS2:` / `w/o MS2:` でない
    # （接頭辞は照合結果ではなく .arf2 の Name に付く）。assign 済みでも元のスペクトルの有無で判断する。
    annotations = load_spot_annotations(str(arf2_path))
    name_by_id = {int(r["MasterAlignmentID"]): r.get("Name") or "" for r in catalog}
    msms = {}
    for sid, ann in annotations.items():
        rep = ann.get("representative") or {}
        msms[sid] = (bool(rep.get("has_msms"))
                     and name_prefix(name_by_id.get(sid, "")) not in ("no MS2", "w/o MS2"))
    return Selection(items=resolved_items, excluded=excluded, groups=resolved_groups,
                     low_reliability=frozenset(low), rows_by_spot=rows_by_spot, msms=msms,
                     catalog_by_id={int(r["MasterAlignmentID"]): r for r in catalog}, caveats=caveats)


def expand_spots(selection: Selection, *, require_msms: bool = False):
    """解決済みの項目をスポット単位に展開する（項目の順 → 表示名の順。同じスポットは 1 回）。"""
    seen = set()
    spots = []
    dropped = []
    for item in selection.items:
        ordered = sorted(item["spots"], key=lambda s: (display_name(s["name"]).casefold(), s["spot_id"]))
        for spot in ordered:
            sid = int(spot["spot_id"])
            if sid in seen:
                continue
            seen.add(sid)
            has_msms = bool(selection.msms.get(sid))
            if require_msms and not has_msms:
                dropped.append(f"#{sid} {spot['name']}")
                continue
            row = selection.catalog_by_id.get(sid, {})
            spots.append({
                "spot_id": sid, "name": spot["name"], "label": display_name(spot["name"]),
                "ontology": spot["ontology"], "item": item["item"],
                "adduct": row.get("AdductType") or "",
                "mz": round(float(row.get("MassCenter") or 0.0), 4),
                "rt": round(float(row.get("RT") or 0.0), 3),
                "msms": has_msms,
            })
    return spots, sorted(dropped)
```

- [ ] **Step 4: 新しいテストが通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_arf_selection.py -q`
Expected: PASS（5 件）

- [ ] **Step 5: `arf_plot_group_intensity` の本体を `build_selection` に置き換える**

`metabolomix/arf/tools.py` の `arf_plot_group_intensity` で、関数内の import 群と `try:` ブロックを次に置き換える（docstring・前提チェック・`arf_state.last_group_intensity = None` は残す）。`from metabolomix.arf import selection as arf_selection` は関数内 import に足す。

```python
    from metabolomix.arf import selection as arf_selection
    from metabolomix.curation import flags as curation_flags
    from metabolomix.msdial import analysis_params
    from metabolomix.plots import group_intensity as gi

    # （前提チェックと last_group_intensity = None は既存のまま）

    try:
        mode = plot_render.resolve_plot_output(output)
        if detection_limit is not None:
            if (isinstance(detection_limit, bool) or not isinstance(detection_limit, (int, float))
                    or not math.isfinite(detection_limit) or detection_limit <= 0):
                raise ValueError(f"detection_limit は 0 より大きい有限の数で指定してください（受け取った値: {detection_limit!r}）。")
        if ncols is not None and (isinstance(ncols, bool) or not isinstance(ncols, int) or ncols < 1):
            raise ValueError(f"ncols は 1 以上の整数で指定してください（受け取った値: {ncols!r}）。")
        selected = arf_selection.build_selection(
            arf_state, arf2_path, items=items, groups=groups,
            low_reliability_samples=low_reliability_samples, apply_curation=apply_curation,
            exclude_auto_likely_wrong=exclude_auto_likely_wrong, standard_samples=standard_samples)
        source = None
        if detection_limit is not None:
            detection_limit, source = float(detection_limit), "argument"
        else:
            param = analysis_params.find_param_file(arf2_path)
            value = analysis_params.read_analysis_params(param)["min_peak_height"] if param else None
            if value:
                detection_limit, source = value, "param_file"
    except curation_flags.FlagFileError as exc:
        return json_payload({"status": "error", "message": str(exc), **exc.details()})
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})

    payload = gi.build_group_intensity_payload(
        selected.items, selected.groups, selected.rows_by_spot, msms=selected.msms,
        low_reliability=selected.low_reliability, excluded=selected.excluded,
        detection_limit=detection_limit, detection_limit_source=source, caveats=selected.caveats)
```

それ以降（`if mode == plot_render.PAYLOAD:` からの戻り値）は既存のまま。不要になった関数内 import（`arf_reader` `arf2_reader` `load_spot_annotations` `name_prefix` `curation_apply` `curation_suggest`）は消す。

- [ ] **Step 6: 群別強度の既存テストが変わらず通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity_tools.py tests/test_group_intensity.py tests/test_arf_selection.py -q`
Expected: PASS（既存テストは無修正）

- [ ] **Step 7: workflow 文書の連鎖を直す**

`docs/workflow/arf.md` の `## arf_plot_group_intensity` の連鎖（手順 5〜18）を次に置き換え、本文の手順番号の言及（「手順 2・8・10・16」「手順 12」「手順 13〜14」「手順 15」「手順 16」「手順 17〜18」「手順 19〜20」）を新しい番号に合わせて書き直す。

```
1. metabolomix/arf/tools.py  arf_plot_group_intensity()
2. └─ metabolomix/plots/render.py  resolve_plot_output()
3. └─ metabolomix/core/mcp_errors.py  missing_state()
4. └─ metabolomix/arf/tools.py  _sibling_arf2_path()
5. └─ metabolomix/arf/selection.py  build_selection()
6. │  └─ metabolomix/arf/reader.py  alignment_feature_row()
7. │  └─ metabolomix/msdial/sample_factors.py  build_sample_facets()
8. │  └─ metabolomix/plots/group_intensity.py  resolve_groups()
9. │  └─ metabolomix/plots/group_intensity.py  resolve_sample_specs()
10.│  └─ metabolomix/arf2/reader.py  load_catalog()
11.│  └─ [apply_curation] metabolomix/curation/apply.py  flags_for_arf2()
12.│  └─ [exclude_auto_likely_wrong] metabolomix/curation/suggest.py  latest_review()
13.│  └─ [standard_samples] metabolomix/plots/group_intensity.py  standard_only_spots()
14.│  └─ metabolomix/plots/group_intensity.py  resolve_items()
15.│  └─ metabolomix/arf2/match_results.py  load_spot_annotations()
16.└─ metabolomix/msdial/analysis_params.py  find_param_file()
17.└─ metabolomix/msdial/analysis_params.py  read_analysis_params()
18.└─ metabolomix/plots/group_intensity.py  build_group_intensity_payload()
19.├─ [output=payload] metabolomix/core/serialization.py  json_payload()
20.└─ [output=image] metabolomix/plots/group_intensity.py  render_group_intensity_plot()
21.   └─ metabolomix/plots/render.py  figure_to_png()
22.   └─ metabolomix/arf/tools.py  _group_intensity_caption()
```

Run: `C:/Python314/python.exe -m pytest tests/test_workflow_docs.py -q`
Expected: PASS

- [ ] **Step 8: コミット（バックグラウンド）**

```bash
git add metabolomix/arf/selection.py metabolomix/arf/tools.py tests/test_arf_selection.py docs/workflow/arf.md
git commit -m "refactor(arf): 項目の解決を arf/selection.py に切り出し、群別強度をそれに載せ替える" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t1.log 2>&1
```

ログ末尾に `passed` と `FAILED` が無いことを確かめる。

---

### Task 2: `save_figure` を作り、旧 4 ツールを統合してから消す

**Files:**
- Modify: `metabolomix/tools/reports.py`
- Modify（文言）: `metabolomix/arf/tools.py` `metabolomix/eic/tools.py` `metabolomix/core/tool_helpers.py` `metabolomix/core/session_state.py` `metabolomix/mztab/dataset_state.py` `metabolomix/plots/volcano.py` `metabolomix/plots/render.py`
- Modify（テスト移行）: `tests/test_report_tools.py` `tests/test_differential_tools.py` `tests/test_eic_plot.py` `tests/test_eic_multi_plot.py` `tests/test_group_intensity_tools.py` `tests/test_missing_state_contract.py`
- Create: `tests/test_save_figure.py`
- 台帳 7 か所（`docs/workflow/plots.md` は全面書き換え）、`docs/output_format/arf.md` `docs/output_format/eic.md` `docs/workflow/arf.md` `docs/workflow/eic.md`

**Interfaces:**
- Produces: `reports.save_figure(kind: str, analysis_id: str, title: str | None = None, source: str = "auto", result_id: str | None = None) -> str`。`reports.FIGURE_KINDS: tuple[str, ...]`、`reports._SOURCE_KINDS`。後続の Task 4 / 10 が `FIGURE_KINDS` に `"species"` / `"pca_loadings"` を足し、`save_figure` に分岐を足す。
- 内部ヘルパ（旧ツールの本体を移す）: `_save_eic(analysis_id, title) -> str`、`_save_group_intensity(analysis_id, title) -> str`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_save_figure.py`:

```python
"""save_figure: 図の保存の一本化。旧 4 ツールと同じ出力になることを縛る。"""
from __future__ import annotations

import json

import pytest

import server
from metabolomix.core import session_state

_PCA_PLOT = {"title": "PCA", "x_label": "PC1 (60.00%)", "y_label": "PC2 (20.00%)",
             "points": [{"x": 1.0, "y": 2.0, "label": "a"}, {"x": -1.0, "y": 0.5, "label": "b"},
                        {"x": 0.2, "y": -1.0, "label": "c"}]}


@pytest.fixture(autouse=True)
def _fresh(tmp_path, monkeypatch):
    session_state.session = server.AnalysisSession()
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    yield


def test_unknown_kind_is_an_error():
    out = json.loads(server.save_figure("heatmap", "x"))
    assert out["status"] == "error" and "kind" in out["message"]


@pytest.mark.parametrize("kind", ["eic", "group_intensity"])
def test_source_is_rejected_for_kinds_without_sources(kind):
    out = json.loads(server.save_figure(kind, "x", source="arf"))
    assert out["status"] == "error" and "source" in out["message"]
    out = json.loads(server.save_figure(kind, "x", result_id="r1"))
    assert out["status"] == "error"


@pytest.mark.parametrize("kind,state,tools", [
    ("pca", "pca_result", ["arf_parser", "arf_pca_preprocessed", "load_dataset", "dataset_pca"]),
    ("volcano", "differential_result", ["arf_differential", "dataset_differential"]),
    ("eic", "eic_plot", ["eic_plot_chromatograms", "eic_plot_compounds"]),
    ("group_intensity", "group_intensity_plot", ["arf_plot_group_intensity"]),
])
def test_missing_state_per_kind(kind, state, tools):
    err = json.loads(server.save_figure(kind, "x"))["error"]
    assert err["code"] == "missing_state" and err["state"] == state and err["required_tools"] == tools


def test_pca_writes_png_and_reports_absolute_path(tmp_path):
    session_state.session.arf.last_pca_plot = dict(_PCA_PLOT)
    msg = server.save_figure("pca", "a-1")
    png = tmp_path / "reports" / "figures" / "a-1_pca.png"
    assert png.is_file() and str(png) in msg and "figures/a-1_pca.png" in msg


def test_group_intensity_writes_png_and_svg(tmp_path):
    payload = {"plot_schema": "lipidmix.group_intensity.v1", "detection_limit": None,
               "groups": [{"label": "g", "samples": ["s1"]}],
               "items": [{"item": "PG", "n_spots": 1, "n_spots_msms": None, "detected": True,
                          "groups": [{"label": "g", "samples": [
                              {"sample": "s1", "value": 100.0, "gap_filled_fraction": 0.0,
                               "low_reliability": False}],
                              "log10_mean": 2.0, "log10_sd": None, "n_in_stats": 1}]}]}
    session_state.session.arf.last_group_intensity = {"payload": payload, "title": None, "ncols": None}
    msg = server.save_figure("group_intensity", "EV membrane")
    pngs = list((tmp_path / "reports" / "figures").glob("*_group_intensity.png"))
    assert len(pngs) == 1 and pngs[0].with_suffix(".svg").is_file()
    assert "group_intensity.png" in msg
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_save_figure.py -q`
Expected: FAIL（`AttributeError: module 'server' has no attribute 'save_figure'`）

- [ ] **Step 3: `reports.py` に `save_figure` を足す（旧ツールはまだ残す）**

`metabolomix/tools/reports.py` の import に `from metabolomix.core.serialization import json_payload` を足し、`__all__` に `"save_figure"` を足す。旧 `save_eic_figure` の本体を `_save_eic(analysis_id, title)`、旧 `save_group_intensity_figure` の本体を `_save_group_intensity(analysis_id, title)` という素の関数（デコレータなし）に移し、旧ツールはそれを呼ぶだけにする。そのうえで次を足す:

```python
#: save_figure が受け付ける図の種類。後続の作業で species / pca_loadings を足す。
FIGURE_KINDS = ("pca", "volcano", "eic", "group_intensity")
#: ARF と mzTab のどちらの結果を描くかを選ぶ必要がある種類（source / result_id を使う）。
_SOURCE_KINDS = ("pca", "volcano")


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
def save_figure(kind: str, analysis_id: str, title: str | None = None,
                source: str = "auto", result_id: str | None = None) -> str:
    """明示的なユーザー要求時だけ、直近の図を reports/figures/<analysis_id>_<kind>.png に保存する。

    kind: "pca"（arf_parser / arf_pca_preprocessed / dataset_pca のスコア）/ "volcano"
      （arf_differential / dataset_differential）/ "eic"（eic_plot_chromatograms /
      eic_plot_compounds）/ "group_intensity"（arf_plot_group_intensity。dpi 300 の PNG と同名 .svg）。
    source / result_id: kind が pca / volcano のときだけ使う。"auto"（既定）はどちらも優先せず、
      有効な結果が 2 つ以上あると AMBIGUOUS_RESULT_SOURCE で止まる（"arf" / "mztab" で指定する）。
    通常の対話描画では呼ばない（各描画ツールが画像を返す）。返り値の相対パスは write_report の
    本文に `![...](figures/<analysis_id>_<kind>.png)` として埋め込める。
    """
    if kind not in FIGURE_KINDS:
        return json_payload({"status": "error", "message":
                             f"kind は {', '.join(FIGURE_KINDS)} のどれかを指定してください（受け取った値: {kind!r}）。"})
    if kind not in _SOURCE_KINDS and (source != "auto" or result_id is not None):
        return json_payload({"status": "error", "message":
                             f"source / result_id は kind が {' / '.join(_SOURCE_KINDS)} のときだけ指定できます"
                             f"（kind={kind!r}）。"})
    if kind == "pca":
        chosen = _choose_figure_result("pca", source, result_id)
        if isinstance(chosen, str):
            return chosen
        return _save_figure(chosen, analysis_id, title, kind="pca", label="PCA")
    if kind == "volcano":
        chosen = _choose_figure_result("differential", source, result_id)
        if isinstance(chosen, str):
            return chosen
        return _save_figure(chosen, analysis_id, title, kind="volcano", label="volcano")
    if kind == "eic":
        return _save_eic(analysis_id, title)
    return _save_group_intensity(analysis_id, title)
```

`server.py` は `from metabolomix.tools.reports import *` なので追加の import は要らない。

- [ ] **Step 4: 新旧が同じ出力になることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_save_figure.py tests/test_report_tools.py tests/test_group_intensity_tools.py tests/test_eic_plot.py tests/test_eic_multi_plot.py tests/test_differential_tools.py -q`
Expected: `tests/test_save_figure.py` が PASS。既存テストも PASS（旧ツールはまだある）。ただし `test_server_registration.py` 等は次の Step で直すまで実行しない。

- [ ] **Step 5: 旧 4 ツールを消す**

`metabolomix/tools/reports.py` から `save_pca_figure` `save_volcano_figure` `save_eic_figure` `save_group_intensity_figure` の関数を削除し、`__all__` から外す。モジュール docstring 1 行目の列挙を `write/read/list_reports, save_figure。` に直す。

- [ ] **Step 6: 既存テストを `save_figure` に移す**

呼び出しだけを置き換える（`def test_save_..._figure` というテスト関数名には当てない）:

```bash
sed -i -E 's/([^_A-Za-z])save_pca_figure\(/\1save_figure("pca", /g; s/([^_A-Za-z])save_volcano_figure\(/\1save_figure("volcano", /g; s/([^_A-Za-z])save_eic_figure\(/\1save_figure("eic", /g; s/([^_A-Za-z])save_group_intensity_figure\(/\1save_figure("group_intensity", /g; s/import save_group_intensity_figure/import save_figure/g' tests/test_report_tools.py tests/test_differential_tools.py tests/test_eic_plot.py tests/test_eic_multi_plot.py tests/test_group_intensity_tools.py tests/test_missing_state_contract.py
grep -n "save_pca_figure\|save_volcano_figure\|save_eic_figure\|save_group_intensity_figure" tests/*.py
```

`grep` に残るのはテスト関数名とコメントだけのはず。`tests/test_report_tools.py` の PCA 案内文の検査（`self.assertIn("save_pca_figure", block)`）は Step 7 の文言変更に合わせて `self.assertIn("save_figure", block)` に直す。`tests/test_missing_state_contract.py` のテスト関数名は `test_save_figure_pca` / `test_save_figure_volcano` / `test_save_figure_group_intensity` / `test_save_figure_eic_offers_both_producers` に改名する。

- [ ] **Step 7: 本体の文言を直す**

旧ツール名を案内している箇所を `save_figure(kind="...")` に直す（意味は変えない）:

- `metabolomix/arf/tools.py`: volcano の案内（`save_volcano_figure を実行します` → `save_figure(kind="volcano") を実行します`、docstring の ``save_volcano_figure`` → ``save_figure(kind="volcano")``、コメント）と群別強度の docstring（`save_group_intensity_figure` → `save_figure(kind="group_intensity")`）・コメント
- `metabolomix/core/tool_helpers.py`: `PNG が必要なときのみ save_pca_figure を実行します。` → `PNG が必要なときのみ save_figure(kind="pca") を実行します。`、docstring の言及
- `metabolomix/eic/tools.py`: ``save_eic_figure`` → ``save_figure(kind="eic")``
- `metabolomix/core/session_state.py`: コメント 3 行（`save_figure が参照`）
- `metabolomix/mztab/dataset_state.py`・`metabolomix/plots/volcano.py`・`metabolomix/plots/render.py`: コメント・docstring の言及

Run: `grep -rn "save_pca_figure\|save_volcano_figure\|save_eic_figure\|save_group_intensity_figure\|save_\*_figure" metabolomix server.py`
Expected: 出力なし

- [ ] **Step 8: 台帳 7 か所を直す（74 → 71）**

1. `tests/test_server_registration.py`: `EXPECTED_TOOLS` から旧 4 つを消し `"save_figure"` を足す。コメントに `# 図の保存の一本化: save_figure に統合し旧 4 ツールを撤去(74→71)。` を足す。
2. `tests/test_tool_annotations.py`: 旧 4 行を消し `"save_figure": LOCAL_WRITE,` を足す。
3. `tests/test_workflow_docs.py`: `"plots.md": ("save_figure",),`、末尾の合計を `50` に。
4. `docs/workflow/plots.md` を書き換える（下記）。
5. `docs/workflow/index.md`: plots.md の行を `| [plots.md](plots.md) | 図の保存（save_figure）と payload 契約の比較 | 1 |`、`合計 50 ツール`、`登録ツール総数は 71（50 + 21）`。
6. `USAGE.md`: 旧 4 行を消して次の 1 行を足し、見出しを `(全71ツール)`、冒頭の流れ図の `save_*_figure` を `save_figure`、本文の `save_eic_figure` の案内を `save_figure(kind="eic")` に直す。

   `| `save_figure` | ユーザーがファイルとしての図を明示的に希望した場合だけ、直近の図を `reports/figures/<analysis_id>_<kind>.png` に保存(`kind`, `analysis_id`, `title`, `source`, `result_id`)。`kind` は `pca`(arf_parser / arf_pca_preprocessed / dataset_pca のスコア)/ `volcano`(arf_differential / dataset_differential。間引き前の全点)/ `eic`(eic_plot_chromatograms / eic_plot_compounds)/ `group_intensity`(arf_plot_group_intensity。dpi 300 の PNG と同名 `.svg`)。`source`/`result_id` は `pca`/`volcano` だけ: 入力元が ARF と mzTab-M の 2 つあり**どちらも優先しない** — 有効な結果が 2 つ以上あると `AMBIGUOUS_RESULT_SOURCE` で止まるので `source`(`arf`/`mztab`)か `result_id` を指定する。前処理をやり直して古くなった結果は選べない。探索専用データセットから描いた図には図の中に但し書きが入る。未知の `kind` と、`source`/`result_id` を使わない `kind` への指定は `{"status":"error"}`。図が無ければ `missing_state`(`required_tools` は kind ごと)。戻り値は保存した絶対パスと埋め込み用の相対パス。 |`
7. `CLAUDE.md`: `ツール 71・リソース 5・リソーステンプレート 3`。

`docs/workflow/plots.md` の新しい内容（payload 契約の表と保存先の説明は旧文書から引き継ぐ）:

```markdown
# ワークフロー: 図の保存と payload 契約

通常の対話では PNG ファイルを作らない。ユーザーが明示的に保存を求めたときだけ `save_figure`
を呼ぶ。`kind` で「どの直近結果を読むか」を振り分け、`session` に記録済みの結果から描き直して
`reports/figures/<analysis_id>_<kind>.png` に書く。

## payload 契約の比較

| payload | 生成元 | 保存 | 特徴 |
|---|---|---|---|
| PCA スコア（`session.arf.last_pca_plot` / `session.dataset.last_pca`） | `arf_parser` / `arf_pca_preprocessed` / `dataset_pca` | `save_figure(kind="pca")` | 散布図 |
| `lipidmix.volcano.v1` | `arf_plot_volcano(output="payload")` | `save_figure(kind="volcano")` | `up`/`down` は全点保持、`ns` のみ間引く。保存は間引き前の全点 |
| `lipidmix.group_intensity.v1` | `arf_plot_group_intensity(output="payload")` | `save_figure(kind="group_intensity")` | 項目 × 群 × 試料の PeakHeight 合計。PNG（dpi 300）に加えて SVG も書く |
| `lipidmix.eic.v1` / `.multi.v1` | `eic_plot_chromatograms` / `eic_plot_compounds(output="payload")` | `save_figure(kind="eic")` | 線グラフ |

描画ツールの**既定は payload ではなく画像**（`output="image"`）。座標点列を LLM の文脈へ流すと
実測で数万トークンかかるのに対し、サーバ側で描いた PNG は 300〜1,000 画像トークンに収まる。
Plotly で描くクライアント（Use-LLLM）は起動 env に `LIPIDMIX_PLOT_OUTPUT=payload` を置くか、
呼び出しごとに `output="payload"` を渡す。判定は `plots/render.py resolve_plot_output()`。
描画は `plots/` 側に置き、画像返しと保存が同じ関数を通る（画面の図とレポートの図がずれない）。

保存先は `_resolve_report_dir()` が解決する。`LIPIDMIX_REPORTS_DIR` の明示先を優先し、書き込めなければ
解析フォルダ配下の `reports/` に退避する。未指定時は解析フォルダの `reports/` を優先し、不可なら
リポジトリ内の `reports/` に退避する。

## save_figure

前提: `kind` に対応する描画・解析を実行済み（無ければ `MissingState`。`required_tools` は kind ごと）
状態変更: `reports/figures/` に PNG（`group_intensity` は SVG も）を書く。

`kind` が `pca` / `volcano` のときだけ入力元を選ぶ（手順 2〜4）。ARF 経路と mzTab-M 経路の結果を
並べ、**どちらも優先しない**——有効な結果が 2 つ以上あれば `AMBIGUOUS_RESULT_SOURCE` で止まり、
`source` か `result_id` の指定を求める。前処理をやり直して古くなった結果は候補にならない。
未知の `kind` と、`source` / `result_id` を使わない `kind` への指定は `{"status": "error"}`。
戻り値には保存した絶対パスを書く（Use-LLLM はこれを拾って画像として表示する）。

1. metabolomix/tools/reports.py  save_figure()
2. ├─ [kind=pca] metabolomix/tools/reports.py  _pca_candidates()
3. ├─ [kind=volcano] metabolomix/tools/reports.py  _differential_candidates()
4. ├─ [kind=pca / volcano] metabolomix/tools/reports.py  _choose_figure_result()
5. │  └─ metabolomix/plots/result_output.py  select_result()
6. │  └─ metabolomix/tools/reports.py  _save_figure()
7. │     └─ metabolomix/plots/result_output.py  save_result_figure()
8. ├─ [kind=eic] metabolomix/tools/reports.py  _save_eic()
9. │  └─ metabolomix/plots/eic.py  render_eic_plot()
10.└─ [kind=group_intensity] metabolomix/tools/reports.py  _save_group_intensity()
11.   └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
```

`docs/output_format/arf.md` §11.5・`docs/output_format/eic.md`・`docs/workflow/arf.md`・`docs/workflow/eic.md` の旧ツール名の言及を `save_figure(kind="...")` に直す。

- [ ] **Step 9: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS

- [ ] **Step 10: コミット（バックグラウンド）**

```bash
git add -A metabolomix tests docs/workflow docs/output_format USAGE.md CLAUDE.md
git commit -m "feat(reports): 図の保存を save_figure に一本化し、旧 4 ツール（save_pca/volcano/eic/group_intensity_figure）を撤去" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t2.log 2>&1
```

---

### Task 3: 分子種ごとの図の payload と描画（`plots/species.py`）

**Files:**
- Create: `metabolomix/plots/species.py`
- Test: `tests/test_species_plot.py`（新規）

**Interfaces:**
- Consumes: `selection.expand_spots` の要素の形（Task 1）
- Produces:
  - `species.SPECIES_SCHEMA = "lipidmix.species_intensity.v1"`、`species.MAX_PANELS = 40`、`species.VALUES = ("share", "height")`
  - `species.build_species_payload(spots, groups, rows_by_spot, *, value, basis_spots=None, basis_items=None, low_reliability=frozenset(), excluded=None, caveats=None) -> dict`
  - `species.render_species_plot(payload, *, title=None, ncols=None) -> matplotlib.figure.Figure`
  - payload のキー: `plot_schema` `value` `stat` `share_basis`（share のとき `{"items", "spot_ids"}`、height のとき `None`）`groups` `low_reliability_samples` `zero_denominator_samples` `spots`（各スポットに `groups`: `[{"label", "samples": [{"sample", "height", "share", "gap_filled", "low_reliability"}], "mean", "sd", "n_in_stats"}]`）`excluded` `caveats`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_species_plot.py`:

```python
"""plots/species.py: 分子種ごとの図の payload と描画（純関数）。"""
from __future__ import annotations

import math

import pytest

from metabolomix.plots import species


def _spot(sid, label, adduct="[M-H]-"):
    return {"spot_id": sid, "name": label, "label": label, "ontology": label.split()[0], "item": label.split()[0],
            "adduct": adduct, "mz": 700.0, "rt": 3.0, "msms": True}


def _rows(values, gap=()):
    return [{"file_name": n, "height": h, "is_gap_filled": n in gap} for n, h in values.items()]


GROUPS = [{"label": "C", "samples": ["c1", "c2"]}, {"label": "P", "samples": ["p1", "p2"]}]
ROWS = {
    0: _rows({"c1": 30, "c2": 20, "p1": 10, "p2": 0}, gap=("p2",)),
    1: _rows({"c1": 70, "c2": 80, "p1": 90, "p2": 0}),
}


def test_share_is_height_over_basis_total_and_zero_total_is_listed():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1"), _spot(1, "PG 16:0_19:1")], GROUPS, ROWS,
                                      value="share")
    c = p["spots"][0]["groups"][0]
    assert [s["share"] for s in c["samples"]] == [30.0, 20.0]
    assert c["mean"] == 25.0 and c["sd"] == pytest.approx(7.0711, abs=1e-4)
    assert p["zero_denominator_samples"] == ["p2"]
    p2 = p["spots"][0]["groups"][1]["samples"][1]
    assert p2["share"] is None and p2["gap_filled"] is True
    assert p["share_basis"]["spot_ids"] == [0, 1]


def test_basis_spots_change_the_denominator():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="share",
                                      basis_spots=[_spot(1, "PG 16:0_19:1")], basis_items=["PG 16:0_19:1"])
    assert p["spots"][0]["groups"][0]["samples"][0]["share"] == pytest.approx(30 / 70 * 100, abs=1e-3)
    assert p["share_basis"] == {"items": ["PG 16:0_19:1"], "spot_ids": [1]}


def test_height_stats_are_log10_and_skip_zero_and_low_reliability():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="height",
                                      low_reliability=frozenset({"c2"}))
    c, pg = p["spots"][0]["groups"]
    assert c["n_in_stats"] == 1 and c["mean"] == pytest.approx(math.log10(30), abs=1e-4) and c["sd"] is None
    assert pg["n_in_stats"] == 1        # p2 は 0 なので統計に入らない
    assert p["share_basis"] is None and p["stat"].startswith("mean ± SD of log10")


def test_group_with_only_low_reliability_samples_has_no_mean():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="share",
                                      low_reliability=frozenset({"p1", "p2"}))
    assert p["spots"][0]["groups"][1]["mean"] is None
    species.render_species_plot(p)      # 落ちない


def test_unknown_value_is_rejected():
    with pytest.raises(ValueError):
        species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="area")


@pytest.mark.parametrize("value", ["share", "height"])
def test_render_draws_one_panel_per_spot_with_adduct_in_title(value):
    import matplotlib.pyplot as plt
    spots = [_spot(0, "DG 16:0_18:1", "[M+NH4]+"), _spot(1, "DG 16:0_18:1", "[M+Na]+")]
    p = species.build_species_payload(spots, GROUPS, ROWS, value=value)
    fig = species.render_species_plot(p, ncols=2)
    titles = [ax.get_title(loc="left") for ax in fig.axes if ax.get_visible() and ax.get_title(loc="left")]
    assert any("[M+NH4]+" in t for t in titles) and any("[M+Na]+" in t for t in titles)
    plt.close(fig)


def test_render_rejects_empty_payload():
    p = species.build_species_payload([], GROUPS, ROWS, value="share")
    with pytest.raises(ValueError):
        species.render_species_plot(p)
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_species_plot.py -q`
Expected: FAIL（`ImportError: cannot import name 'species'`）

- [ ] **Step 3: `metabolomix/plots/species.py` を書く**

```python
"""分子種ごとの図: 1 パネル = 1 アラインメントスポット、群ごとに試料別の値を並べる（MCP 非依存）。

value="share": 試料ごとに「スポットの高さ ÷ 分母（share_basis のスポットの高さの合計）× 100」。
  棒 = 群平均、ひげ = SD（n ≥ 2）、点 = 各試料。
value="height": PeakHeight。log10 の平均 ± SD（群別強度と同じ規則）、0 は軸の下端。
縦軸はパネルごと（分子種ごとに桁が違い、群間の変動を見るのが目的のため）。検定はしない。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §4。
"""
from __future__ import annotations

import math
import statistics

from metabolomix.plots.group_intensity import GROUP_COLORS

SPECIES_SCHEMA = "lipidmix.species_intensity.v1"
MAX_PANELS = 40
VALUES = ("share", "height")


def build_species_payload(spots, groups, rows_by_spot, *, value, basis_spots=None, basis_items=None,
                          low_reliability=frozenset(), excluded=None, caveats=None):
    """スポット × 群 × 試料の値を組み立てる（描画しない）。"""
    if value not in VALUES:
        raise ValueError(f"value は {' / '.join(VALUES)} のどちらかを指定してください（受け取った値: {value!r}）。")
    index = {sid: {r.get("file_name"): r for r in rows} for sid, rows in rows_by_spot.items()}
    samples = []
    for group in groups:
        for name in group["samples"]:
            if name not in samples:
                samples.append(name)

    def height(sid, name):
        row = index.get(sid, {}).get(name)
        return float(row.get("height") or 0.0) if row else 0.0

    denominators = {}
    basis_ids = []
    zero = []
    if value == "share":
        basis_ids = [int(s["spot_id"]) for s in (basis_spots if basis_spots is not None else spots)]
        denominators = {name: sum(height(sid, name) for sid in basis_ids) for name in samples}
        zero = [name for name in samples if denominators[name] <= 0]

    out = []
    for spot in spots:
        sid = int(spot["spot_id"])
        out_groups = []
        for group in groups:
            entries = []
            for name in group["samples"]:
                row = index.get(sid, {}).get(name)
                h = height(sid, name)
                share = None
                if value == "share" and denominators.get(name, 0.0) > 0:
                    share = round(h / denominators[name] * 100.0, 4)
                entries.append({"sample": name, "height": round(h, 1), "share": share,
                                "gap_filled": bool(row.get("is_gap_filled")) if row else False,
                                "low_reliability": name in low_reliability})
            if value == "share":
                stats = [e["share"] for e in entries if e["share"] is not None and not e["low_reliability"]]
            else:
                stats = [math.log10(e["height"]) for e in entries if e["height"] > 0 and not e["low_reliability"]]
            out_groups.append({
                "label": group["label"], "samples": entries,
                "mean": round(statistics.mean(stats), 4) if stats else None,
                "sd": round(statistics.stdev(stats), 4) if len(stats) > 1 else None,
                "n_in_stats": len(stats),
            })
        out.append({**spot, "groups": out_groups})
    return {
        "plot_schema": SPECIES_SCHEMA,
        "value": value,
        "stat": "mean ± SD of share (%)" if value == "share" else "mean ± SD of log10(peak height)",
        "share_basis": ({"items": list(basis_items) if basis_items else None, "spot_ids": basis_ids}
                        if value == "share" else None),
        "groups": [{"label": g["label"], "samples": list(g["samples"])} for g in groups],
        "low_reliability_samples": sorted(low_reliability),
        "zero_denominator_samples": zero,
        "spots": out,
        "excluded": {key: list(values) for key, values in (excluded or {}).items()},
        "caveats": list(caveats or []),
    }


def render_species_plot(payload, *, title=None, ncols=None):
    """payload を 1 スポット 1 パネルの図にする。"""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    spots = payload["spots"]
    if not spots:
        raise ValueError("描く分子種がありません。")
    labels = [g["label"] for g in payload["groups"]]
    share = payload["value"] == "share"
    ncols = ncols or min(5, len(spots))
    nrows = math.ceil(len(spots) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.3 * ncols, 2.25 * nrows + 0.7), squeeze=False)
    for ax, spot in zip(axes.flat, spots):
        has_zero = False
        for i, group in enumerate(spot["groups"]):
            color = GROUP_COLORS[i % len(GROUP_COLORS)]
            if group["mean"] is not None:
                if share:
                    ax.bar(i, group["mean"], width=0.62, color=color, alpha=0.55, edgecolor=color, zorder=2)
                else:
                    ax.hlines(group["mean"], i - 0.28, i + 0.28, color="black", linewidth=1.6, zorder=3)
                if group["sd"] is not None:
                    ax.errorbar(i, group["mean"], yerr=group["sd"], color="black", capsize=3, linewidth=0.9, zorder=3)
            n = len(group["samples"])
            for j, s in enumerate(group["samples"]):
                x = i + (j - (n - 1) / 2) * 0.12
                if share:
                    if s["share"] is None:
                        continue
                    y = s["share"]
                elif s["height"] <= 0:
                    has_zero = True
                    ax.scatter(x, 0.0, marker="v", s=18, facecolor="white", edgecolor="#555555", zorder=4)
                    continue
                else:
                    y = math.log10(s["height"])
                ax.scatter(x, y, s=18, marker="D" if s["gap_filled"] else "o",
                           facecolor="white" if s["low_reliability"] else color,
                           edgecolor="black", linewidth=0.5, zorder=4)
        ax.set_title(f"{spot['label']}\n{spot['adduct']} (ID {spot['spot_id']})", fontsize=7, loc="left")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=6)
        ax.set_xlim(-0.6, len(labels) - 0.4)
        ax.tick_params(axis="y", labelsize=6)
        if share:
            ax.set_ylim(bottom=0)
        elif has_zero:
            ax.set_ylim(bottom=-0.3)
        ax.set_ylabel("% of basis" if share else "log10(peak height)", fontsize=6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in list(axes.flat)[len(spots):]:
        ax.axis("off")
    legend = [Patch(facecolor=GROUP_COLORS[i % len(GROUP_COLORS)], alpha=0.7, label=label)
              for i, label in enumerate(labels)]
    legend += [
        Line2D([], [], marker="o", ls="", mfc="#888888", mec="black", label="detected peak"),
        Line2D([], [], marker="D", ls="", mfc="#888888", mec="black", label="gap-filled value"),
        Line2D([], [], marker="o", ls="", mfc="white", mec="black", label="low-reliability (not in mean/SD)"),
    ]
    if not share:
        legend.append(Line2D([], [], marker="v", ls="", mfc="white", mec="#555555", label="zero (axis floor)"))
    fig.legend(handles=legend, loc="lower center", ncol=min(len(legend), 6), fontsize=6.5, frameon=False)
    default = ("Share of each species (bar = mean, whisker = SD)" if share
               else "Peak height of each species (line = mean of log10, whisker = SD)")
    fig.suptitle(title or default, fontsize=10)
    fig.tight_layout(rect=(0, 0.6 / fig.get_figheight(), 1, 0.97))
    return fig
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_species_plot.py -q`
Expected: PASS（8 件）

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/plots/species.py tests/test_species_plot.py
git commit -m "feat(plots): 分子種ごとの図（割合 % / PeakHeight）の payload と描画" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t3.log 2>&1
```

---

### Task 4: `arf_plot_species` ツールと `save_figure(kind="species")`

**Files:**
- Create: `metabolomix/arf/species_tools.py`
- Modify: `server.py`（import 1 行）、`metabolomix/core/session_state.py`（`ArfState.last_species_plot`）、`metabolomix/tools/reports.py`（kind `species`）
- Test: `tests/test_species_tools.py`（新規）、`tests/test_save_figure.py`、`tests/test_missing_state_contract.py`
- 台帳 7 か所（`docs/workflow/arf.md` に節）、`docs/output_format/arf.md`（payload の節）

**Interfaces:**
- Consumes: `selection.build_selection` / `expand_spots`（Task 1）、`species.build_species_payload` / `render_species_plot` / `MAX_PANELS`（Task 3）、`reports.FIGURE_KINDS`（Task 2）
- Produces: `species_tools.arf_plot_species(...)`、`species_tools._require_arf_with_arf2() -> Path | str`（Task 7 も使う。`str` は missing_state の封筒）、`session.arf.last_species_plot = {"payload", "title", "ncols"}`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_species_tools.py`（Task 7 でこのファイルに PCA のテストを足す）:

```python
"""arf_plot_species / arf_pca_species（ツール層）。fixture はテストが作る。"""
from __future__ import annotations

import json

import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state
from tests.curation_fixtures import arf2_spot_raw, arf_row, match_result, write_arf, write_arf2

STEM = "AlignmentResult_2026_01_01_00_00_00"
SAMPLES = ["blank_1", "ctrl_1", "ctrl_2", "ctrl_3", "ko_1", "ko_2", "ko_3"]


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
        _group([10, 5000, 6000, 5500, 2000, 2500, 2200]),     # 0 PG 16:0_18:1（MS/MS）
        _group([10, 3000, 3500, 3200, 4000, 4200, 4100]),     # 1 PG 16:0_19:1（MS/MS なし）
        _group([0, 2000, 2100, 1900, 4000, 4300, 300]),       # 2 DGDG 16:0_18:1（MS/MS）
        _group([0, 900, 950, 920, 800, 700, 760]),            # 3 DG 16:0_18:1 [M+NH4]+（MS/MS）
        _group([0, 300, 320, 310, 280, 260, 270]),            # 4 DG 16:0_18:1 [M+Na]+（MS/MS）
    ])
    arf2 = write_arf2(tmp_path / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=0, name="PG 34:1|PG 16:0_18:1", ontology="PG", adduct="[M-H]-",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=1, name="PG 35:1|PG 16:0_19:1", ontology="PG", adduct="[M-H]-"),
        arf2_spot_raw(spot_id=2, name="DGDG 34:1|DGDG 16:0_18:1", ontology="DGDG", adduct="[M+CH3COO]-",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=3, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+NH4]+",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=4, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+Na]+",
                      matches=[match_result()]),
    ])
    session_state.session = session_state.AnalysisSession()
    with open(arf, "rb") as handle:
        session_state.session.arf.features = arf_reader.deserialize(handle)
    session_state.session.arf.current_file_path = str(arf)
    session_state.session.arf.class_index = None
    return {"arf": arf, "arf2": arf2}


def _species(**kw):
    from metabolomix.arf.species_tools import arf_plot_species
    kw.setdefault("groups", ["ctrl", "ko"])
    return json.loads(arf_plot_species(output="payload", **kw))


def test_missing_arf_returns_missing_state():
    from metabolomix.arf.species_tools import arf_plot_species
    session_state.session = session_state.AnalysisSession()
    out = json.loads(arf_plot_species(items=["PG"], groups=["ctrl"], output="payload"))
    assert out["error"]["code"] == "missing_state"
    assert out["error"]["required_tools"] == ["arf_parser", "load_dataset"]


def test_share_payload_expands_classes_and_uses_items_as_default_basis(loaded):
    p = _species(items=["PG", "DGDG"])
    assert p["plot_schema"] == "lipidmix.species_intensity.v1" and p["value"] == "share"
    assert [s["spot_id"] for s in p["spots"]] == [0, 1, 2]
    assert p["share_basis"]["spot_ids"] == [0, 1, 2]
    c1 = p["spots"][0]["groups"][0]["samples"][0]
    assert c1["share"] == pytest.approx(5000 / (5000 + 3000 + 2000) * 100, abs=1e-3)


def test_share_basis_and_require_msms(loaded):
    p = _species(items=["PG"], share_basis=["PG", "DGDG"], require_msms=True)
    assert [s["spot_id"] for s in p["spots"]] == [0]
    assert p["share_basis"]["spot_ids"] == [0, 2]        # 分母にも MS/MS の条件が効く
    assert p["excluded"]["no_msms"] == ["#1 PG 35:1|PG 16:0_19:1"]


def test_basis_without_the_plotted_species_adds_caveat(loaded):
    p = _species(items=["PG"], share_basis=["DGDG"])
    assert any("分母に含まれない" in c for c in p["caveats"])


def test_same_name_different_adduct_are_separate_panels(loaded):
    p = _species(items=["DG"], value="height")
    assert [(s["spot_id"], s["adduct"]) for s in p["spots"]] == [(3, "[M+NH4]+"), (4, "[M+Na]+")]


def test_too_many_panels_is_an_error_and_clears_slot(loaded, monkeypatch):
    from metabolomix.plots import species
    monkeypatch.setattr(species, "MAX_PANELS", 2)
    _species(items=["PG"])
    assert session_state.session.arf.last_species_plot is not None
    out = _species(items=["PG", "DGDG"])
    assert out["status"] == "error" and "40" not in out["message"] and "2" in out["message"]
    assert session_state.session.arf.last_species_plot is None


@pytest.mark.parametrize("kw", [{"value": "area"}, {"ncols": 0}, {"groups": ["nonexistent"]}])
def test_invalid_input_returns_error(loaded, kw):
    out = _species(items=["PG"], **kw)
    assert out["status"] == "error" and out["message"]


def test_image_output_caption_and_session(loaded):
    from mcp.server.fastmcp import Image
    from metabolomix.arf.species_tools import arf_plot_species
    out = arf_plot_species(items=["PG"], groups=["ctrl", "ko"], low_reliability_samples=["ko_3"])
    assert isinstance(out, list) and isinstance(out[1], Image)
    assert "2 パネル" in out[0] and "ctrl n=3" in out[0]
    assert session_state.session.arf.last_species_plot["payload"]["spots"][0]["spot_id"] == 0


def test_save_figure_species_writes_png_and_svg(loaded, tmp_path, monkeypatch):
    import server
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    _species(items=["PG"])
    msg = server.save_figure("species", "EV-PG")
    png = tmp_path / "reports" / "figures" / "EV-PG_species.png"
    assert png.is_file() and png.with_suffix(".svg").is_file() and str(png) in msg
```

`tests/test_save_figure.py` の `test_missing_state_per_kind` の parametrize に `("species", "species_plot", ["arf_plot_species"]),` を足し、`test_source_is_rejected_for_kinds_without_sources` の parametrize に `"species"` を足す。`tests/test_missing_state_contract.py` に次を足す:

```python
    def test_arf_plot_species(self):
        self.assert_missing(
            server.arf_plot_species(items=["PG"], groups=["x"], output="payload"),
            "arf_dataset",
            ["arf_parser", "load_dataset"],
            "load_dataset",
        )
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_species_tools.py tests/test_save_figure.py -q`
Expected: FAIL（`ModuleNotFoundError: No module named 'metabolomix.arf.species_tools'`）

- [ ] **Step 3: セッションの欄を足す**

`metabolomix/core/session_state.py` の `ArfState.__init__` の「直近解析の描画用データ」に `self.last_species_plot = None  # arf_plot_species / save_figure(kind="species") が参照` と `self.last_species_pca = None  # arf_pca_species / save_figure(kind="pca", source="species") / plot_pca_loadings が参照` を足し、`reset_analysis()` にも同じ 2 つを `None` にする行を足す。

- [ ] **Step 4: `metabolomix/arf/species_tools.py` を書く（この段階では `arf_plot_species` だけ）**

```python
"""ARF の分子種ごとの図（arf_plot_species）と分子種 PCA（arf_pca_species）の MCP ツール。

項目の解決は arf/selection.py、描画は plots/species.py / plots/pca_scores.py。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §4・§5。
"""
from __future__ import annotations

from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.arf import selection as arf_selection
from metabolomix.core import mcp_errors, session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.curation import flags as curation_flags
from metabolomix.plots import render as plot_render
from metabolomix.plots import species as species_plot

__all__ = ["arf_plot_species"]

_EXCLUDED_LABELS = (("internal_standard", "内部標準"), ("manual", "手動除外"), ("curation", "判断"),
                    ("auto_likely_wrong", "自動判定"), ("standard_only", "標準液"), ("no_msms", "MS/MS なし"))


def _require_arf_with_arf2():
    """ARF と兄弟 .arf2 があればその .arf2 のパス、無ければ missing_state の封筒（str）。"""
    from metabolomix.arf.tools import _sibling_arf2_path

    arf_state = session_state.session.arf
    if arf_state.features is None or not str(arf_state.current_file_path or "").lower().endswith(".arf"):
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"], "先に load_dataset で ARF データを読み込んでください。")
    arf2_path = _sibling_arf2_path()
    if arf2_path is None:
        return mcp_errors.missing_state(
            "sibling_arf2", ["arf_parser", "load_dataset"],
            "同じアラインメントの .arf2 が見つかりません（名前とクラスの解決に要ります）。")
    return arf2_path


def _check_ncols(ncols):
    if ncols is not None and (isinstance(ncols, bool) or not isinstance(ncols, int) or ncols < 1):
        raise ValueError(f"ncols は 1 以上の整数で指定してください（受け取った値: {ncols!r}）。")


def _excluded_summary(excluded: dict) -> str:
    return "・".join(f"{label} {len(excluded.get(key, []))}" for key, label in _EXCLUDED_LABELS)


def _species_caption(payload: dict) -> str:
    groups = "、".join(f"{g['label']} n={len(g['samples'])}" for g in payload["groups"])
    kind = "割合 %" if payload["value"] == "share" else "PeakHeight"
    caption = (f"分子種ごとの{kind}: {len(payload['spots'])} パネル（{groups}）。"
               f"除外 — {_excluded_summary(payload['excluded'])}。")
    if payload["zero_denominator_samples"]:
        caption += f" 分母が 0 で描かなかった試料: {'、'.join(payload['zero_denominator_samples'])}。"
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return caption


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_plot_species(
    items: list[str],
    groups: list[str],
    value: str = "share",
    share_basis: list[str] | None = None,
    low_reliability_samples: list[str] | None = None,
    apply_curation: bool = True,
    exclude_auto_likely_wrong: bool = False,
    standard_samples: list[str] | None = None,
    require_msms: bool = False,
    title: str | None = None,
    ncols: int | None = None,
    output: str | None = None,
) -> list | str:
    """選んだクラスを分子種（アラインメントスポット）に展開し、1 分子種 1 パネルで群ごとの試料別の値を並べる。

    - items: クラス名（.arf2 の Ontology に完全一致）は配下のスポットに展開、それ以外は分子種名。
      同じ名前でも付加イオンが違えば別パネル（見出しに付加イオンと ID）。1 回 40 パネルまで。
    - groups: arf_plot_group_intensity と同じトークン規則（QC・ブランクは qc / blank を書いたときだけ）。
    - value: "share"（既定。試料ごとに 高さ ÷ 分母の合計 × 100。棒 = 群平均、ひげ = SD）/ "height"
      （PeakHeight を log10 で。横線 = log10 の平均、ひげ = SD、0 は軸の下端）。検定はしない。
    - share_basis: 割合の分母にする項目（省略時は items）。例: PG のページでも分母は
      ["PG", "DGDG", "MGDG", "CL"]。分母にも同じ除外を適用する。
    - 除外は arf_plot_group_intensity と同じ（内部標準・standard_samples だけにある分子種・
      apply_curation の判断・exclude_auto_likely_wrong・arf_exclude したスポット）。require_msms=True で
      MS/MS の裏付けがない分子種も除く。除いたものは payload の excluded。
    - low_reliability_samples: 白抜きで描き、平均と SD から外す。菱形は gap-filled の値。
    - output: "image"（既定）/ "payload"（lipidmix.species_intensity.v1）。PNG ファイルは
      save_figure(kind="species")。
    - 入力の誤りは {"status": "error", "message": ...}。失敗した呼び出しの後は前回の図を破棄する。
    """
    arf2_path = _require_arf_with_arf2()
    if isinstance(arf2_path, str):
        return arf2_path
    arf_state = session_state.session.arf
    arf_state.last_species_plot = None
    filters = dict(low_reliability_samples=low_reliability_samples, apply_curation=apply_curation,
                   exclude_auto_likely_wrong=exclude_auto_likely_wrong, standard_samples=standard_samples)
    try:
        mode = plot_render.resolve_plot_output(output)
        if value not in species_plot.VALUES:
            raise ValueError(f"value は share / height のどちらかを指定してください（受け取った値: {value!r}）。")
        _check_ncols(ncols)
        selected = arf_selection.build_selection(arf_state, arf2_path, items=items, groups=groups, **filters)
        spots, no_msms = arf_selection.expand_spots(selected, require_msms=require_msms)
        if not spots:
            raise ValueError("描く分子種がありません（項目が当たらないか、すべて除外されました）。")
        if len(spots) > species_plot.MAX_PANELS:
            raise ValueError(f"パネルが {len(spots)} 枚になります（1 回 {species_plot.MAX_PANELS} 枚まで）。"
                             "items を分けて呼んでください（share_basis を同じにすれば割合はそろいます）。")
        basis_spots = None
        caveats = list(selected.caveats)
        if value == "share" and share_basis:
            basis_sel = arf_selection.build_selection(arf_state, arf2_path, items=share_basis, groups=groups, **filters)
            basis_spots, _ = arf_selection.expand_spots(basis_sel, require_msms=require_msms)
            outside = sorted({s["spot_id"] for s in spots} - {s["spot_id"] for s in basis_spots})
            if outside:
                caveats.append(f"分母（share_basis）に含まれない分子種が {len(outside)} 個あります（ID {outside}）。"
                               "割合が 100 % を超えることがあります。")
    except curation_flags.FlagFileError as exc:
        return json_payload({"status": "error", "message": str(exc), **exc.details()})
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})

    excluded = {**selected.excluded, "no_msms": no_msms}
    payload = species_plot.build_species_payload(
        spots, selected.groups, selected.rows_by_spot, value=value, basis_spots=basis_spots,
        basis_items=share_basis if value == "share" else None, low_reliability=selected.low_reliability,
        excluded=excluded, caveats=caveats)
    if mode == plot_render.PAYLOAD:
        arf_state.last_species_plot = {"payload": payload, "title": title, "ncols": ncols}
        return json_payload(payload)
    png = plot_render.figure_to_png(species_plot.render_species_plot(payload, title=title, ncols=ncols))
    arf_state.last_species_plot = {"payload": payload, "title": title, "ncols": ncols}
    return [_species_caption(payload), Image(data=png, format="png")]
```

`server.py` の `from metabolomix.arf.tools import *` の次の行に `from metabolomix.arf.species_tools import *  # arf_plot_species, arf_pca_species` を足す。

- [ ] **Step 5: `save_figure` に `species` を足す**

`metabolomix/tools/reports.py`:

```python
FIGURE_KINDS = ("pca", "volcano", "eic", "group_intensity", "species")
```

`save_figure` の docstring の kind の列挙に `"species"（arf_plot_species。dpi 300 の PNG と同名 .svg）` を足し、末尾の `return _save_group_intensity(...)` を次に置き換える:

```python
    if kind == "group_intensity":
        return _save_group_intensity(analysis_id, title)
    return _save_species(analysis_id, title)


def _save_species(analysis_id: str, title: str | None) -> str:
    from metabolomix.plots.species import render_species_plot

    last = getattr(session_state.session.arf, "last_species_plot", None)
    if not last or not last.get("payload"):
        return mcp_errors.missing_state(
            "species_plot", ["arf_plot_species"],
            "先に arf_plot_species を実行してください（分子種ごとの図がありません）。")
    out_path = _figures_dir() / f"{knowledge_store.make_slug(analysis_id)}_species.png"
    fig = render_species_plot(last["payload"], title=title or last.get("title"), ncols=last.get("ncols"))
    _write_png_and_svg(fig, out_path)
    return (f"分子種ごとの図を保存: {out_path}（同名の .svg も保存）\n"
            f"本文に ![species](figures/{out_path.name}) で埋め込めます。")
```

`_save_group_intensity` と共通の 2 つのヘルパを足し、`_save_group_intensity` もこれを使うように直す:

```python
def _figures_dir():
    figures_dir = _resolve_report_dir() / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    return figures_dir


def _write_png_and_svg(fig, out_path) -> None:
    """報告書用の図を dpi 300 の PNG と同名 SVG で書く（Figure は閉じる）。"""
    try:
        fig.savefig(out_path, dpi=300, format="png", bbox_inches="tight")
        fig.savefig(out_path.with_suffix(".svg"), format="svg", bbox_inches="tight")
    finally:
        plt.close(fig)
```

- [ ] **Step 6: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_species_tools.py tests/test_save_figure.py tests/test_group_intensity_tools.py -q`
Expected: PASS

- [ ] **Step 7: 台帳 7 か所を直す（71 → 72）**

1. `EXPECTED_TOOLS` に `"arf_plot_species"`、コメント `# 分子種ごとの図: arf_plot_species を追加(71→72)。`
2. `tests/test_tool_annotations.py`: `"arf_plot_species": READ_ONLY,`
3. `tests/test_workflow_docs.py`: `"arf.md"` のタプルに `"arf_plot_species"`、合計 `51`
4. `docs/workflow/arf.md` の末尾に節を足す:

```markdown
## arf_plot_species

前提: ARF を読み込み済みで、同じアラインメントの `.arf2` が隣にある（無ければ `MissingState`）
状態変更: `session.arf.last_species_plot`（payload・title・ncols）を更新。ファイルは書かない。
`save_figure(kind="species")`（[plots.md](plots.md)）がここから読んで保存する。

選んだクラスを分子種（アラインメントスポット）に展開し、1 分子種 1 パネルで群ごとの試料別の値を並べる。
項目・群・除外の解決は `arf_plot_group_intensity` と同じ `arf/selection.py`。`value="share"` の分母は
`share_basis`（省略時は `items`）を同じ規則で解決したスポットの合計（手順 6 をもう 1 回呼ぶ）。
入力の誤り・パネル超過（40 枚）・判断の記録ファイルの破損は `{"status":"error","message":...}` を返して終わる。

1. metabolomix/arf/species_tools.py  arf_plot_species()
2. └─ metabolomix/arf/species_tools.py  _require_arf_with_arf2()
3. │  └─ metabolomix/core/mcp_errors.py  missing_state()
4. │  └─ metabolomix/arf/tools.py  _sibling_arf2_path()
5. └─ metabolomix/plots/render.py  resolve_plot_output()
6. └─ metabolomix/arf/selection.py  build_selection()
7. │  └─ metabolomix/plots/group_intensity.py  resolve_items()
8. └─ metabolomix/arf/selection.py  expand_spots()
9. └─ metabolomix/plots/species.py  build_species_payload()
10.├─ [output=payload] metabolomix/core/serialization.py  json_payload()
11.└─ [output=image] metabolomix/plots/species.py  render_species_plot()
12.   └─ metabolomix/plots/render.py  figure_to_png()
13.   └─ metabolomix/arf/species_tools.py  _species_caption()
```

   `docs/workflow/plots.md` の payload 表に行 `| lipidmix.species_intensity.v1 | arf_plot_species(output="payload") | save_figure(kind="species") | スポット × 群 × 試料の高さと割合。PNG（dpi 300）と SVG |` を足し、`## save_figure` の連鎖の末尾を次にする:

```
10.├─ [kind=group_intensity] metabolomix/tools/reports.py  _save_group_intensity()
11.│  └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
12.└─ [kind=species] metabolomix/tools/reports.py  _save_species()
13.   └─ metabolomix/plots/species.py  render_species_plot()
```

5. `docs/workflow/index.md`: arf.md の件数 `12`、`合計 51 ツール`、`登録ツール総数は 72（51 + 21）`
6. `USAGE.md`: 見出し `(全72ツール)`、`arf_plot_group_intensity` の行の後に:

   `| `arf_plot_species` | 選んだクラスを分子種(アラインメントスポット)に展開し、1 分子種 1 パネルで群ごとの試料別の値を並べる(`items`, `groups`, `value`, `share_basis`, `low_reliability_samples`, `apply_curation`, `exclude_auto_likely_wrong`, `standard_samples`, `require_msms`, `title`, `ncols`, `output`)。同じ名前でも付加イオンが違えば別パネル、1 回 40 パネルまで。`value="share"`(既定)は試料ごとに 高さ ÷ 分母 × 100(分母は `share_basis`、省略時は `items`)で棒 = 群平均・ひげ = SD、`value="height"` は log10(PeakHeight) の点と平均 ± SD。除外は `arf_plot_group_intensity` と同じで、`require_msms=True` で MS/MS の裏付けがない分子種も除く(`excluded` の `no_msms`)。`low_reliability_samples` は白抜きで平均・SD から外す。既定は PNG 画像と 1 行の説明、`output="payload"` で `lipidmix.species_intensity.v1`。前提: ARF の読み込みと兄弟 `.arf2`(無ければ `missing_state`)。入力の誤りは `{"status":"error","message":...}`。ファイルは `save_figure(kind="species")`。 |`

   `save_figure` の行の kind の列挙に `/ `species`(arf_plot_species。dpi 300 の PNG と同名 `.svg`)` を足す。
7. `CLAUDE.md`: `ツール 72・リソース 5・リソーステンプレート 3`

`docs/output_format/arf.md` の §11.5 の後に `### 11.6 `arf_plot_species` の返り値（`lipidmix.species_intensity.v1`）` を足し、payload の各キー（Task 3 の Interfaces の一覧）と意味（`share` は分母に対する %、分母が 0 の試料は `share=null` で `zero_denominator_samples` に入る、`mean`/`sd` は `stat` の定義に従う、`excluded.no_msms`）を箇条書きで書く。

- [ ] **Step 8: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS

- [ ] **Step 9: コミット（バックグラウンド）**

```bash
git add -A metabolomix server.py tests docs/workflow docs/output_format USAGE.md CLAUDE.md
git commit -m "feat(arf): 分子種ごとの図 arf_plot_species と save_figure(kind=species) を足す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t4.log 2>&1
```

---

### Task 5: 計算に使う試料を選べる PCA と相関ローディング（`analysis/pca.py`）

**Files:**
- Modify: `metabolomix/analysis/pca.py`
- Test: `tests/test_analysis_pca.py`

**Interfaces:**
- Produces:
  - `pca.run_pca_fit_subset(matrix, fit_mask, n_components=None) -> dict`（`components`（全試料のスコア）、`explained_variance_ratio`、`singular_values`、`loadings`（主成分 × 残した特徴量）、`kept_features`（元の列の番号）、`n_fit`）
  - `pca.loading_correlations(loadings, singular_values, n_fit) -> list[list[float]]`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_analysis_pca.py` の末尾に足す:

```python
def test_run_pca_fit_subset_projects_unfitted_rows_and_drops_constant_columns():
    import numpy as np
    from metabolomix.analysis.pca import run_pca_fit_subset

    rng = np.random.default_rng(0)
    fit_rows = rng.normal(size=(6, 4))
    fit_rows[:, 3] = 5.0                                  # 計算に使う試料で一定 → 外す
    outlier = np.array([[10.0, -10.0, 3.0, 99.0]])
    matrix = np.vstack([fit_rows, outlier])
    fit = np.array([True] * 6 + [False])
    out = run_pca_fit_subset(matrix, fit, n_components=2)
    assert out["kept_features"] == [0, 1, 2] and out["n_fit"] == 6
    assert len(out["components"]) == 7                    # 投影した試料にもスコアがある
    # 計算に使った試料だけで決めた軸: 外れ値を入れても fit 試料のスコアは変わらない
    alone = run_pca_fit_subset(fit_rows, np.ones(6, bool), n_components=2)
    assert np.allclose(np.abs(out["components"][:6]), np.abs(alone["components"]), atol=1e-9)


def test_loading_correlations_equal_pearson_r_on_fitted_rows():
    import numpy as np
    from metabolomix.analysis.pca import loading_correlations, run_pca_fit_subset

    rng = np.random.default_rng(1)
    matrix = rng.normal(size=(8, 5))
    fit = np.array([True] * 7 + [False])
    out = run_pca_fit_subset(matrix, fit, n_components=2)
    r = np.asarray(loading_correlations(out["loadings"], out["singular_values"], out["n_fit"]))
    scores = np.asarray(out["components"])[fit]
    for k in range(2):
        for j in range(5):
            expected = np.corrcoef(matrix[fit][:, j], scores[:, k])[0, 1]
            assert r[k, j] == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("fit", [[True, True, False, False], [True, True, True, False]])
def test_run_pca_fit_subset_needs_three_fitted_rows_and_two_varying_columns(fit):
    import numpy as np
    from metabolomix.analysis.pca import run_pca_fit_subset

    matrix = np.array([[1.0, 2.0], [2.0, 2.0], [3.0, 2.0], [4.0, 2.0]])   # 列 1 は一定
    with pytest.raises(ValueError):
        run_pca_fit_subset(matrix, np.array(fit))
```

（`tests/test_analysis_pca.py` の先頭に `import pytest` が無ければ足す。）

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_analysis_pca.py -q`
Expected: FAIL（`ImportError: cannot import name 'run_pca_fit_subset'`）

- [ ] **Step 3: 実装する**

`metabolomix/analysis/pca.py` の先頭に `import math` を足し、末尾に足す:

```python
def run_pca_fit_subset(matrix, fit_mask, n_components: int | None = None) -> dict:
    """計算に使う試料（fit_mask=True）だけで autoscale と PCA の軸を決め、全試料をその軸へ投影する。

    低信頼の試料を主成分の計算から外しつつ、どこに落ちるかは見たいときに使う。autoscale は
    計算に使う試料の平均と**母標準偏差**（ddof=0。StandardScaler と同じ）。計算に使う試料で
    分散が 0 の列は外し、`kept_features` に残した列の番号を返す。
    """
    matrix = np.asarray(matrix, dtype=float)
    fit = np.asarray(fit_mask, dtype=bool)
    if matrix.ndim != 2 or fit.shape != (matrix.shape[0],):
        raise ValueError("matrix は 2 次元、fit_mask は行数と同じ長さにしてください。")
    n_fit = int(fit.sum())
    if n_fit < 3:
        raise ValueError(f"PCA の計算に使う試料は 3 以上必要です（現在: {n_fit}）。")
    fit_rows = matrix[fit]
    sd = fit_rows.std(axis=0)
    kept = np.flatnonzero(sd > 0)
    if kept.size < 2:
        raise ValueError(f"計算に使う試料で値が変わる特徴量が 2 未満のため PCA を計算できません（現在: {kept.size}）。")
    mean = fit_rows[:, kept].mean(axis=0)
    scaled = (matrix[:, kept] - mean) / sd[kept]
    max_components = min(n_fit, kept.size)
    target = max_components if n_components is None else min(n_components, max_components)
    pca = PCA(n_components=target).fit(scaled[fit])
    return {
        "components": pca.transform(scaled).tolist(),
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "singular_values": pca.singular_values_.tolist(),
        "loadings": pca.components_.tolist(),
        "kept_features": kept.tolist(),
        "n_fit": n_fit,
    }


def loading_correlations(loadings, singular_values, n_fit: int) -> list[list[float]]:
    """各特徴量と主成分スコアの相関 r（主成分 × 特徴量）。

    r_jk = 固有ベクトルの成分_jk × 特異値_k / √n。母標準偏差で autoscale した PCA でだけ
    成り立つ（n は主成分の計算に使った試料数）。中心化だけの PCA には使わない。
    """
    coef = np.asarray(loadings, dtype=float)
    s = np.asarray(singular_values, dtype=float)[: coef.shape[0]]
    return (coef * s[:, None] / math.sqrt(n_fit)).tolist()
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_analysis_pca.py -q`
Expected: PASS

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/analysis/pca.py tests/test_analysis_pca.py
git commit -m "feat(analysis): 計算に使う試料を選べる PCA と相関ローディングを足す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t5.log 2>&1
```

---

### Task 6: PCA スコア図の描画を共通化（`plots/pca_scores.py`）

**Files:**
- Create: `metabolomix/plots/pca_scores.py`
- Modify: `metabolomix/plots/result_output.py`（`_render` の pca 分岐）
- Test: `tests/test_pca_scores.py`（新規）、既存 `tests/test_result_output.py`

**Interfaces:**
- Produces: `pca_scores.render_pca_scores(plot: dict, *, title: str | None = None) -> Figure`。`plot` は `{"title", "x_label", "y_label", "points": [{"x", "y", "label", "group"?, "fitted"?}]}`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_pca_scores.py`:

```python
"""plots/pca_scores.py: PCA スコア図（群の色分け・投影点の白抜き・枠外の矢印）。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from metabolomix.plots.pca_scores import render_pca_scores


def _plot(points):
    return {"title": "t", "x_label": "PC1 (50.0%)", "y_label": "PC2 (20.0%)", "points": points}


def test_groups_get_a_legend_and_projected_points_are_open():
    fig = render_pca_scores(_plot([
        {"x": 1, "y": 1, "label": "c1", "group": "C"}, {"x": 1.2, "y": 0.8, "label": "c2", "group": "C"},
        {"x": -1, "y": -1, "label": "p1", "group": "P"}, {"x": -1.1, "y": -0.9, "label": "p2", "group": "P"},
        {"x": -0.5, "y": 0.2, "label": "p3", "group": "P", "fitted": False},
    ]))
    texts = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert "C" in texts and "P" in texts and any("projected" in t for t in texts)
    plt.close(fig)


def test_projected_point_outside_the_fitted_range_is_drawn_as_an_arrow():
    fig = render_pca_scores(_plot([
        {"x": 1, "y": 1, "label": "a", "group": "C"}, {"x": -1, "y": -1, "label": "b", "group": "C"},
        {"x": 0.5, "y": -0.5, "label": "c", "group": "P"}, {"x": 40, "y": -60, "label": "far", "group": "P", "fitted": False},
    ]))
    notes = [t.get_text() for t in fig.axes[0].texts]
    assert any("far" in n and "off-scale" in n for n in notes)
    plt.close(fig)


def test_points_without_groups_still_render():
    fig = render_pca_scores(_plot([{"x": 1, "y": 2, "label": "a"}, {"x": -1, "y": 0, "label": "b"}]), title="X")
    assert fig.axes[0].get_title() == "X"
    plt.close(fig)


def test_empty_points_are_rejected():
    with pytest.raises(ValueError):
        render_pca_scores(_plot([]))
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_scores.py -q`
Expected: FAIL（`ModuleNotFoundError`）

- [ ] **Step 3: `metabolomix/plots/pca_scores.py` を書く**

```python
"""PCA スコア図（MCP 非依存）。分子種 PCA の画像返しと save_figure(kind="pca") が共有する。

points の `group` があれば色分けと凡例、`fitted=False`（主成分の計算に使わず投影だけした試料）は
白抜き。軸の範囲は計算に使った点に合わせ、そこから外れた投影点は枠の端に矢印と座標で示す
（外れ値 1 点のために他の点が潰れないように）。
"""
from __future__ import annotations

from metabolomix.plots.group_intensity import GROUP_COLORS


def render_pca_scores(plot: dict, *, title: str | None = None):
    import matplotlib.pyplot as plt

    points = plot.get("points") or []
    if not points:
        raise ValueError("PCA 図に描ける点がありません。")
    fitted = [p for p in points if p.get("fitted", True)] or points
    xs = [float(p["x"]) for p in fitted]
    ys = [float(p["y"]) for p in fitted]
    pad_x = 0.15 * (max(xs) - min(xs)) + 0.5
    pad_y = 0.15 * (max(ys) - min(ys)) + 0.5
    lo_x, hi_x, lo_y, hi_y = min(xs) - pad_x, max(xs) + pad_x, min(ys) - pad_y, max(ys) + pad_y
    groups = []
    for p in points:
        g = p.get("group")
        if g is not None and g not in groups:
            groups.append(g)
    color_of = {g: GROUP_COLORS[(i + 1) % len(GROUP_COLORS)] for i, g in enumerate(groups)}

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(lo_y, hi_y)
    labelled = set()
    has_projected = False
    for p in points:
        x, y = float(p["x"]), float(p["y"])
        projected = not p.get("fitted", True)
        has_projected = has_projected or projected
        if projected and not (lo_x <= x <= hi_x and lo_y <= y <= hi_y):
            cx, cy = min(max(x, lo_x), hi_x), min(max(y, lo_y), hi_y)
            mx, my = (lo_x + hi_x) / 2, (lo_y + hi_y) / 2
            ax.annotate(f"{p.get('label', '')} (off-scale: {x:.0f}, {y:.0f})", xy=(cx, cy),
                        xytext=(cx - 0.25 * (cx - mx), cy - 0.25 * (cy - my)), fontsize=7, ha="center",
                        arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "#555555"})
            continue
        g = p.get("group")
        label = g if (g is not None and g not in labelled and not projected) else None
        if label:
            labelled.add(g)
        ax.scatter(x, y, s=60, facecolor="white" if projected else color_of.get(g, "#3182bd"),
                   edgecolor="black", linewidth=0.7, zorder=3, label=label)
        if p.get("label"):
            ax.annotate(str(p["label"]), (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.axhline(0, color="#eeeeee", lw=0.6, zorder=1)
    ax.axvline(0, color="#eeeeee", lw=0.6, zorder=1)
    ax.set_xlabel(plot.get("x_label", "PC1"))
    ax.set_ylabel(plot.get("y_label", "PC2"))
    ax.set_title(title or plot.get("title", "PCA"))
    if groups or has_projected:
        handles, labels = ax.get_legend_handles_labels()
        if has_projected:
            handles.append(plt.Line2D([], [], marker="o", ls="", mfc="white", mec="black"))
            labels.append("projected (not used to fit)")
        ax.legend(handles, labels, fontsize=7, loc="best")
    fig.tight_layout()
    return fig
```

- [ ] **Step 4: `save_result_figure` の PCA をこれに載せ替える**

`metabolomix/plots/result_output.py` の `_render` の `if kind == "pca":` 分岐を次に置き換え、import を `from metabolomix.core.tool_helpers import dataset_pca_plot` と `from metabolomix.plots.pca_scores import render_pca_scores` に直す（`_pca_scatter_arrays` の import は消す）:

```python
    if kind == "pca":
        # 描画は arf_pca_species（画像返し）と共有する。群があれば色分けされる。
        return render_pca_scores(_pca_plot_payload(result), title=title)
```

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_scores.py tests/test_result_output.py tests/test_report_tools.py tests/test_save_figure.py tests/test_pipeline_end_to_end.py -q`
Expected: PASS

- [ ] **Step 6: workflow 文書を直す**

`docs/workflow/plots.md` の `## save_figure` の連鎖の手順 7 の下に `   └─ metabolomix/plots/pca_scores.py  render_pca_scores()` を足す（`[kind=pca]` のとき）。

- [ ] **Step 7: コミット（バックグラウンド）**

```bash
git add metabolomix/plots/pca_scores.py metabolomix/plots/result_output.py tests/test_pca_scores.py docs/workflow/plots.md
git commit -m "feat(plots): PCA スコア図の描画を pca_scores に共通化（群の色分け・投影点）" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t6.log 2>&1
```

---

### Task 7: `arf_pca_species` ツールと PCA の保存候補 `source="species"`

**Files:**
- Modify: `metabolomix/arf/species_tools.py`、`metabolomix/tools/reports.py`（`_pca_candidates`、`save_figure` の docstring）
- Test: `tests/test_species_tools.py`、`tests/test_missing_state_contract.py`
- 台帳 7 か所（`docs/workflow/arf.md`）、`docs/output_format/arf.md`

**Interfaces:**
- Consumes: `run_pca_fit_subset` / `loading_correlations`（Task 5）、`render_pca_scores`（Task 6）、`_require_arf_with_arf2` / `build_selection` / `expand_spots`
- Produces: `species_tools.arf_pca_species(...)`、`session.arf.last_species_pca` = `{"title", "x_label", "y_label", "points": [{"x", "y", "label", "group", "fitted"}], "explained_variance_ratio", "loadings_rows": [{"feature_id", "label", "ontology", "adduct", "mz", "rt", "coefficient": [..], "r": [..]}], "scaling": "autoscale", "n_fit", "n_species", "dropped_zero_variance", "settings", "provenance"}`。Task 10 は `loadings_rows` と `explained_variance_ratio` を読む。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_species_tools.py` の末尾に足す:

```python
def _pca(**kw):
    from metabolomix.arf.species_tools import arf_pca_species
    kw.setdefault("groups", ["ctrl", "ko"])
    return json.loads(arf_pca_species(output="payload", **kw))


def test_pca_species_projects_low_reliability_and_orients(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], low_reliability_samples=["ko_3"], orient_by="ko")
    pts = {pt["label"]: pt for pt in p["points"]}
    assert pts["20261006_ko_3"]["fitted"] is False and pts["20261006_ctrl_1"]["fitted"] is True
    ko_x = [pts[f"20261006_ko_{i}"]["x"] for i in (1, 2)]
    assert sum(ko_x) / 2 > 0                                    # orient_by の群が正
    last = session_state.session.arf.last_species_pca
    assert last["n_fit"] == 5 and last["scaling"] == "autoscale"
    assert {row["feature_id"] for row in last["loadings_rows"]} == {"0", "1", "2", "3", "4"}
    assert all(-1.0001 <= v <= 1.0001 for row in last["loadings_rows"] for v in row["r"])


def test_pca_species_total_normalization_with_zero_total_names_the_sample(loaded):
    session_state.session.arf.excluded_spots = set()
    out = _pca(items=["DG"], groups=["blank", "ctrl"], normalize="total")
    assert out["status"] == "error" and "blank_1" in out["message"]


def test_pca_species_sample_in_two_groups_is_counted_once(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], groups=["ctrl", "ctrl_1", "ko"])
    labels = [pt["label"] for pt in p["points"]]
    assert labels.count("20261006_ctrl_1") == 1


@pytest.mark.parametrize("kw", [{"normalize": "median"}, {"orient_by": "nope"}, {"items": ["PG 16:0_18:1"]},
                                {"groups": ["ctrl"], "low_reliability_samples": ["ctrl_1"]}])
def test_pca_species_invalid_input_returns_error_and_clears_slot(loaded, kw):
    _pca(items=["PG", "DGDG", "DG"])
    assert session_state.session.arf.last_species_pca is not None
    out = _pca(**{"items": ["PG", "DGDG", "DG"], **kw})
    assert out["status"] == "error" and out["message"]
    assert session_state.session.arf.last_species_pca is None


def test_pca_species_image_and_save_figure_needs_source_when_ambiguous(loaded, tmp_path, monkeypatch):
    from mcp.server.fastmcp import Image
    import server
    from metabolomix.arf.species_tools import arf_pca_species
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    out = arf_pca_species(items=["PG", "DGDG", "DG"], groups=["ctrl", "ko"])
    assert isinstance(out[1], Image) and "PC1" in out[0]
    session_state.session.arf.last_pca_plot = {"title": "x", "x_label": "PC1", "y_label": "PC2",
                                               "points": [{"x": 0, "y": 0, "label": "a"}]}
    ambiguous = json.loads(server.save_figure("pca", "s1"))
    assert ambiguous["error"]["code"] == "AMBIGUOUS_RESULT_SOURCE"
    msg = server.save_figure("pca", "s1", source="species")
    assert (tmp_path / "reports" / "figures" / "s1_pca.png").is_file() and "source=species" in msg
```

`tests/test_missing_state_contract.py` に足す:

```python
    def test_arf_pca_species(self):
        self.assert_missing(
            server.arf_pca_species(items=["PG"], groups=["x"], output="payload"),
            "arf_dataset",
            ["arf_parser", "load_dataset"],
            "load_dataset",
        )
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_species_tools.py -q`
Expected: FAIL（`ImportError: cannot import name 'arf_pca_species'`）

- [ ] **Step 3: `arf_pca_species` を書く**

`metabolomix/arf/species_tools.py` の `__all__` を `["arf_plot_species", "arf_pca_species"]` にし、import に `import numpy as np`、`from metabolomix.analysis.pca import loading_correlations, run_pca_fit_subset`、`from metabolomix.analysis.result_state import array_fingerprint, new_provenance`、`from metabolomix.plots.pca_scores import render_pca_scores` を足して、末尾に:

```python
_PCA_COMPONENTS = 5


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_pca_species(
    items: list[str],
    groups: list[str],
    low_reliability_samples: list[str] | None = None,
    apply_curation: bool = True,
    exclude_auto_likely_wrong: bool = False,
    standard_samples: list[str] | None = None,
    require_msms: bool = False,
    normalize: str = "none",
    log_transform: bool = True,
    orient_by: str | None = None,
    title: str | None = None,
    output: str | None = None,
) -> list | str:
    """選んだ分子種（arf_plot_species と同じ項目指定・除外）で PCA を回し、スコア図を返す。

    - 対象の試料は groups に入る試料（2 群に入る試料は最初の群の色で 1 回だけ）。
    - low_reliability_samples は主成分の計算に使わず、投影だけする（図では白抜き）。
    - normalize: "none"（既定）/ "total"（試料ごとに選んだ分子種の合計で割り、全試料の合計の中央値を掛ける。
      試料量の差を除いて組成を比べるとき）。log_transform（既定 True）で log10(x + 1)。
      その後、計算に使う試料の平均と母標準偏差で autoscale。分散 0 の分子種は外す。
    - orient_by: 群名。その群の平均スコアが正になるよう主成分の符号をそろえる（符号は本来任意）。
    - ローディング（固有ベクトルの成分と相関 r）はセッションに残し、plot_pca_loadings(source="species")
      が描く。スコア図のファイルは save_figure(kind="pca", source="species")。
    - output: "image"（既定）/ "payload"（点列・寄与率・PC1 の r の上位と下位）。
    - 計算に使う試料が 3 未満などの入力の誤りは {"status": "error", "message": ...}。
    """
    arf2_path = _require_arf_with_arf2()
    if isinstance(arf2_path, str):
        return arf2_path
    arf_state = session_state.session.arf
    arf_state.last_species_pca = None
    settings = {"items": list(items), "groups": list(groups), "normalize": normalize,
                "log_transform": log_transform, "require_msms": require_msms, "orient_by": orient_by,
                "low_reliability_samples": list(low_reliability_samples or [])}
    try:
        mode = plot_render.resolve_plot_output(output)
        if normalize not in ("none", "total"):
            raise ValueError(f"normalize は none / total のどちらかを指定してください（受け取った値: {normalize!r}）。")
        selected = arf_selection.build_selection(
            arf_state, arf2_path, items=items, groups=groups, low_reliability_samples=low_reliability_samples,
            apply_curation=apply_curation, exclude_auto_likely_wrong=exclude_auto_likely_wrong,
            standard_samples=standard_samples)
        labels = [g["label"] for g in selected.groups]
        if orient_by is not None and orient_by not in labels:
            raise ValueError(f"orient_by は groups のどれかを指定してください（{labels}）。")
        spots, no_msms = arf_selection.expand_spots(selected, require_msms=require_msms)
        if len(spots) < 2:
            raise ValueError(f"PCA には分子種が 2 つ以上必要です（現在: {len(spots)}）。")
        samples, group_of = [], {}
        for group in selected.groups:
            for name in group["samples"]:
                if name not in group_of:
                    group_of[name] = group["label"]
                    samples.append(name)
        index = {s["spot_id"]: {r.get("file_name"): r for r in selected.rows_by_spot.get(s["spot_id"], [])}
                 for s in spots}
        X = np.array([[float((index[s["spot_id"]].get(n) or {}).get("height") or 0.0) for s in spots]
                      for n in samples])
        if normalize == "total":
            totals = X.sum(axis=1)
            zero = [n for n, t in zip(samples, totals) if t <= 0]
            if zero:
                raise ValueError(f"選んだ分子種の合計が 0 の試料があり、合計で割れません: {'、'.join(zero)}")
            X = X / totals[:, None] * float(np.median(totals))
        if log_transform:
            X = np.log10(X + 1.0)
        fit = np.array([n not in selected.low_reliability for n in samples])
        result = run_pca_fit_subset(X, fit, n_components=_PCA_COMPONENTS)
    except curation_flags.FlagFileError as exc:
        return json_payload({"status": "error", "message": str(exc), **exc.details()})
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})

    scores = np.asarray(result["components"])
    coef = np.asarray(result["loadings"])
    if orient_by is not None:
        mask = np.array([group_of[n] == orient_by for n in samples]) & fit
        for k in range(coef.shape[0]):
            if scores[mask, k].mean() < 0:
                scores[:, k] *= -1
                coef[k] *= -1
    r = np.asarray(loading_correlations(coef, result["singular_values"], result["n_fit"]))
    kept = [spots[i] for i in result["kept_features"]]
    evr = result["explained_variance_ratio"]
    loadings_rows = [{
        "feature_id": str(s["spot_id"]), "label": s["label"], "ontology": s["ontology"], "adduct": s["adduct"],
        "mz": s["mz"], "rt": s["rt"],
        "coefficient": [round(float(coef[k][i]), 4) for k in range(coef.shape[0])],
        "r": [round(float(r[k][i]), 4) for k in range(coef.shape[0])],
    } for i, s in enumerate(kept)]
    plot = {
        "title": title or "PCA of selected species",
        "x_label": f"PC1 ({evr[0] * 100:.1f}%)",
        "y_label": f"PC2 ({evr[1] * 100:.1f}%)" if len(evr) > 1 else "PC2",
        "points": [{"x": round(float(scores[i][0]), 4),
                    "y": round(float(scores[i][1]), 4) if scores.shape[1] > 1 else 0.0,
                    "label": n, "group": group_of[n], "fitted": bool(fit[i])} for i, n in enumerate(samples)],
        "explained_variance_ratio": [round(float(v), 4) for v in evr],
        "loadings_rows": loadings_rows, "scaling": "autoscale", "n_fit": result["n_fit"],
        "n_species": len(kept), "dropped_zero_variance": len(spots) - len(kept),
        "excluded": {**selected.excluded, "no_msms": no_msms}, "caveats": list(selected.caveats),
        "settings": settings,
        "provenance": new_provenance(arf_state, kind="pca", input_fingerprint=array_fingerprint(X),
                                     effective_parameters=settings),
    }
    ranked = sorted(loadings_rows, key=lambda row: row["r"][0], reverse=True)
    top = [f"{row['label']} {row['r'][0]:+.2f}" for row in ranked[:5]]
    bottom = [f"{row['label']} {row['r'][0]:+.2f}" for row in ranked[-5:][::-1]]
    caption = (f"分子種 PCA: {len(kept)} 分子種、計算に使った試料 {result['n_fit']}"
               f"（投影 {int((~fit).sum())}）。PC1 {evr[0] * 100:.1f}%"
               + (f"、PC2 {evr[1] * 100:.1f}%" if len(evr) > 1 else "") + "。"
               f" PC1 の r 上位: {', '.join(top)}。下位: {', '.join(bottom)}。"
               + (f" 分散 0 で外した分子種 {plot['dropped_zero_variance']}。" if plot["dropped_zero_variance"] else "")
               + " ローディング図は plot_pca_loadings(source=\"species\")。")
    if mode == plot_render.PAYLOAD:
        arf_state.last_species_pca = plot
        return json_payload({k: plot[k] for k in ("title", "x_label", "y_label", "points",
                                                  "explained_variance_ratio", "n_species", "n_fit",
                                                  "dropped_zero_variance", "excluded", "caveats")}
                            | {"pc1_r_top": top, "pc1_r_bottom": bottom,
                               "result_id": plot["provenance"].get("result_id")})
    png = plot_render.figure_to_png(render_pca_scores(plot, title=plot["title"]))
    arf_state.last_species_pca = plot
    return [caption, Image(data=png, format="png")]
```

- [ ] **Step 4: PCA の保存候補に `species` を足す**

`metabolomix/tools/reports.py` の `_pca_candidates()` の ARF の候補の後に:

```python
    species = getattr(session_state.session.arf, "last_species_pca", None)
    if species and species.get("points"):
        prov = species.get("provenance") or {}
        candidates.append({
            "source": "species", "result_id": prov.get("result_id"),
            "dataset_id": prov.get("dataset_id"), "valid": True,
            "result": species, "ds": None})
```

`_figure_missing_state("pca")` の `required_tools` に `"arf_pca_species"` を `"arf_pca_preprocessed"` の後へ足し、`save_figure` の docstring の pca の説明に `arf_pca_species` と `source="species"` を足す。`tests/test_save_figure.py` の pca の期待値も `["arf_parser", "arf_pca_preprocessed", "arf_pca_species", "load_dataset", "dataset_pca"]` に直す（`tests/test_missing_state_contract.py` の `test_save_figure_pca` も同じ）。

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_species_tools.py tests/test_save_figure.py tests/test_missing_state_contract.py -q`
Expected: PASS

- [ ] **Step 6: 台帳 7 か所を直す（72 → 73）**

1. `EXPECTED_TOOLS` に `"arf_pca_species"`、コメント `# 分子種 PCA: arf_pca_species を追加(72→73)。`
2. `"arf_pca_species": READ_ONLY,`
3. `"arf.md"` のタプルに `"arf_pca_species"`、合計 `52`
4. `docs/workflow/arf.md` の末尾に節:

```markdown
## arf_pca_species

前提: ARF を読み込み済みで、同じアラインメントの `.arf2` が隣にある（無ければ `MissingState`）
状態変更: `session.arf.last_species_pca`（スコア・ローディング全量・寄与率・provenance）を更新。ファイルは書かない。
`save_figure(kind="pca", source="species")` と `plot_pca_loadings(source="species")` がここから読む。

分子種の選び方は `arf_plot_species` と同じ（手順 4〜5）。試料 × 分子種の PeakHeight 行列を作り、
`normalize="total"` なら試料ごとの合計で割り、log10(x + 1)、計算に使う試料（低信頼を除く）で
autoscale して PCA（手順 6）。低信頼の試料は投影だけする。`orient_by` で符号をそろえ、相関 r を計算する（手順 7）。
入力の誤り（計算に使う試料が 3 未満・合計 0 の試料・未知の `orient_by` など）は `{"status":"error"}`。

1. metabolomix/arf/species_tools.py  arf_pca_species()
2. └─ metabolomix/arf/species_tools.py  _require_arf_with_arf2()
3. └─ metabolomix/plots/render.py  resolve_plot_output()
4. └─ metabolomix/arf/selection.py  build_selection()
5. └─ metabolomix/arf/selection.py  expand_spots()
6. └─ metabolomix/analysis/pca.py  run_pca_fit_subset()
7. └─ metabolomix/analysis/pca.py  loading_correlations()
8. └─ metabolomix/analysis/result_state.py  new_provenance()
9. ├─ [output=payload] metabolomix/core/serialization.py  json_payload()
10.└─ [output=image] metabolomix/plots/pca_scores.py  render_pca_scores()
11.   └─ metabolomix/plots/render.py  figure_to_png()
```

   `docs/workflow/plots.md` の payload 表の PCA スコアの行の生成元に `arf_pca_species`（`session.arf.last_species_pca`）を足す。
5. `docs/workflow/index.md`: arf.md の件数 `13`、`合計 52 ツール`、`登録ツール総数は 73（52 + 21）`
6. `USAGE.md`: 見出し `(全73ツール)`、`arf_plot_species` の行の後に:

   `| `arf_pca_species` | 選んだ分子種(`arf_plot_species` と同じ項目指定・除外)で PCA を回し、スコア図を返す(`items`, `groups`, `low_reliability_samples`, `apply_curation`, `exclude_auto_likely_wrong`, `standard_samples`, `require_msms`, `normalize`, `log_transform`, `orient_by`, `title`, `output`)。対象は `groups` の試料。`low_reliability_samples` は主成分の計算に使わず投影だけする(白抜き)。`normalize="total"` は試料ごとに選んだ分子種の合計で割る(試料量の差を除く)。`log_transform`(既定 True)で log10(x+1)、計算に使う試料で autoscale、分散 0 の分子種は外す。`orient_by` の群の平均スコアが正になるよう符号をそろえる。ローディング(成分と相関 r)はセッションに残し `plot_pca_loadings(source="species")` が描く。スコア図のファイルは `save_figure(kind="pca", source="species")`。既定は PNG 画像と説明(寄与率・PC1 の r の上位と下位)、`output="payload"` で点列と要約。前提: ARF と兄弟 `.arf2`(無ければ `missing_state`)。計算に使う試料が 3 未満・合計 0 の試料などは `{"status":"error","message":...}`。 |`

   `save_figure` の行の `source` の説明に `species`(arf_pca_species)を足す。
7. `CLAUDE.md`: `ツール 73・リソース 5・リソーステンプレート 3`

`docs/output_format/arf.md` に `### 11.7 `arf_pca_species` の返り値` を足す（点列の `fitted`、寄与率、`pc1_r_top` / `pc1_r_bottom` の意味、r は計算に使った試料での相関で −1〜1、符号は `orient_by` で決めない限り任意、低信頼の試料のスコアは投影）。

- [ ] **Step 7: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS

- [ ] **Step 8: コミット（バックグラウンド）**

```bash
git add -A metabolomix tests docs/workflow docs/output_format USAGE.md CLAUDE.md
git commit -m "feat(arf): 分子種を選んで回す PCA arf_pca_species と PCA 保存候補 source=species を足す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t7.log 2>&1
```

---

### Task 8: 既存の PCA がローディングを残すようにする（ARF・mzTab）

**Files:**
- Modify: `metabolomix/core/tool_helpers.py`（`_remember_arf_pca_plot`）、`metabolomix/arf/tools.py`（2 か所の呼び出し）、`metabolomix/analysis/dataset_analysis.py`（`run_dataset_pca`）、`metabolomix/tools/dataset_analysis_tools.py`（要約から外すキー）
- Test: `tests/test_report_tools.py` か新規 `tests/test_pca_loadings_storage.py`

**Interfaces:**
- Produces:
  - `session.arf.last_pca_plot` に追加されるキー: `loadings`（主成分 × 特徴量）、`singular_values`、`explained_variance_ratio`、`feature_names`（`Spot_<id>_<prop>`）、`n_fit`（試料数）、`scaling`（`"autoscale"`）
  - `ds.last_pca` に追加されるキー: `singular_values`（`dataset_pca` の戻り値の要約には載せない）

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_pca_loadings_storage.py`:

```python
"""既存の PCA（ARF / mzTab）がローディング図に要る情報をセッションに残すこと。"""
from __future__ import annotations

import numpy as np

from metabolomix.core import session_state
from metabolomix.core.tool_helpers import _remember_arf_pca_plot


def test_remember_arf_pca_plot_keeps_loadings_and_feature_names():
    session_state.session = session_state.AnalysisSession()
    pca = {"components": [[1.0, 0.5], [-1.0, -0.5], [0.2, 0.1]], "explained_variance_ratio": [0.7, 0.2],
           "singular_values": [3.0, 1.0], "loadings": [[0.6, 0.8], [0.8, -0.6]]}
    _remember_arf_pca_plot(pca, ["a", "b", "c"], title="t", feature_names=["Spot_1_height", "Spot_2_height"])
    plot = session_state.session.arf.last_pca_plot
    assert plot["loadings"] == [[0.6, 0.8], [0.8, -0.6]] and plot["singular_values"] == [3.0, 1.0]
    assert plot["feature_names"] == ["Spot_1_height", "Spot_2_height"]
    assert plot["n_fit"] == 3 and plot["scaling"] == "autoscale"
    assert plot["explained_variance_ratio"] == [0.7, 0.2]
    assert len(plot["points"]) == 3                        # 既存のスコアはそのまま


def test_remember_arf_pca_plot_without_feature_names_stores_none():
    session_state.session = session_state.AnalysisSession()
    pca = {"components": [[1.0, 0.5], [-1.0, -0.5]], "explained_variance_ratio": [0.7, 0.3],
           "singular_values": [3.0, 1.0], "loadings": [[1.0, 0.0], [0.0, 1.0]]}
    _remember_arf_pca_plot(pca, ["a", "b"], title="t")
    assert session_state.session.arf.last_pca_plot["feature_names"] is None


def test_run_dataset_pca_returns_singular_values():
    from metabolomix.analysis.dataset_analysis import run_dataset_pca

    class DS:
        pp_matrix = np.array([[1.0, 2.0, 3.0], [2.0, 1.0, 0.5], [3.0, 5.0, 1.0], [0.5, 0.2, 2.0]])
        pp_sample_names = ["s1", "s2", "s3", "s4"]
        pp_feature_names = ["f1", "f2", "f3"]
        roles = {}
        sample_meta = {}

    out = run_dataset_pca(DS(), n_components=2)
    assert len(out["singular_values"]) == 2 and len(out["loadings"]) == 2
```

（`run_dataset_pca` が読むのは `pp_matrix` `pp_sample_names` `roles` `sample_meta` だけ。`_require_pp_matrix` は `pp_matrix` が None なら止まる。）

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_loadings_storage.py -q`
Expected: FAIL（`TypeError: _remember_arf_pca_plot() got an unexpected keyword argument 'feature_names'`）

- [ ] **Step 3: 実装する**

`metabolomix/core/tool_helpers.py` の `_remember_arf_pca_plot` に引数 `feature_names: list[str] | None = None` を足し、`arf.last_pca_plot = {...}` の dict に次を足す（provenance の前）:

```python
        # ローディング図（plot_pca_loadings）用。LLM への戻り値には載せない。
        # r（相関）は scaling が autoscale のときだけ 成分 × 特異値 / √n で出せる（run_pca は既定 autoscale）。
        "loadings": pca_result.get("loadings"),
        "singular_values": pca_result.get("singular_values"),
        "explained_variance_ratio": pca_result.get("explained_variance_ratio"),
        "feature_names": list(feature_names) if feature_names is not None else None,
        "n_fit": len(sample_names),
        "scaling": "autoscale",
```

`metabolomix/arf/tools.py` の 2 か所の `_remember_arf_pca_plot(...)` 呼び出しに `feature_names=feature_names` を足す（`arf_parser` と `arf_pca_preprocessed`。どちらも同じスコープに `feature_names` がある）。

`metabolomix/analysis/dataset_analysis.py` の `run_dataset_pca` の戻り値 dict に `"singular_values": pca["singular_values"],` を `"loadings"` の隣に足す。`metabolomix/tools/dataset_analysis_tools.py` の要約から外すキーの集合 `("loadings", "provenance")` を `("loadings", "singular_values", "provenance")` にする。

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_loadings_storage.py tests/test_report_tools.py tests/test_dataset_analysis_tools.py -q`
（`tests/test_dataset_analysis_tools.py` が無ければ `tests -q -k "pca"`）
Expected: PASS

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/core/tool_helpers.py metabolomix/arf/tools.py metabolomix/analysis/dataset_analysis.py metabolomix/tools/dataset_analysis_tools.py tests/test_pca_loadings_storage.py
git commit -m "feat(pca): ARF と mzTab の PCA がローディング図に要る情報をセッションに残す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t8.log 2>&1
```

---

### Task 9: ローディング図の組み立てと描画（`plots/pca_loadings.py`）

**Files:**
- Create: `metabolomix/plots/pca_loadings.py`
- Test: `tests/test_pca_loadings_plot.py`（新規）

**Interfaces:**
- Produces:
  - `pca_loadings.LOADINGS_SCHEMA = "lipidmix.pca_loadings.v1"`、`pca_loadings.MAX_ALL_FEATURES = 60`
  - 特徴量の共通の形 `feature = {"feature_id": str, "label": str, "ontology": str | None, "adduct": str, "mz": float | None, "rt": float | None, "coefficient": list[float], "r": list[float] | None}`（`coefficient` / `r` は主成分 1, 2, … の順）
  - `pca_loadings.build_loadings_payload(features, *, source, result_id, explained_variance_ratio, pcs, top_n, value, r_available, caveats=None) -> dict`
  - `pca_loadings.render_loadings_plot(payload, *, title=None) -> Figure`
  - payload のキー: `plot_schema` `source` `result_id` `value`（実際に描いた値の種類）`value_requested` `top_n` `aligned` `panels`（`[{"pc": "PC1", "explained_pct", "rows": [{"feature_id", "label", "ontology", "adduct", "mz", "rt", "coefficient", "r", "value"}]}]`）`caveats`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_pca_loadings_plot.py`:

```python
"""plots/pca_loadings.py: ローディングの選び方と描画（純関数）。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from metabolomix.plots import pca_loadings as pl


def _feat(i, c1, c2, ontology="PG"):
    return {"feature_id": str(i), "label": f"PG {i}:0", "ontology": ontology, "adduct": "[M-H]-",
            "mz": 700.0, "rt": 3.0, "coefficient": [c1, c2], "r": [c1 * 2, c2 * 2]}


FEATS = [_feat(i, c1, c2) for i, (c1, c2) in enumerate(
    [(0.5, 0.1), (0.4, -0.2), (0.1, 0.3), (-0.1, 0.05), (-0.3, -0.4), (-0.45, 0.2)])]


def _payload(**kw):
    args = dict(source="species", result_id="r1", explained_variance_ratio=[0.6, 0.2], pcs=[1, 2],
                top_n=2, value="r", r_available=True)
    args.update(kw)
    return pl.build_loadings_payload(FEATS, **args)


def test_top_n_takes_positive_and_negative_per_pc():
    p = _payload()
    pc1 = [row["feature_id"] for row in p["panels"][0]["rows"]]
    assert pc1 == ["0", "1", "4", "5"] or set(pc1) == {"0", "1", "4", "5"}
    assert p["panels"][0]["explained_pct"] == 60.0 and p["aligned"] is False
    assert p["panels"][0]["rows"][0]["value"] == p["panels"][0]["rows"][0]["r"]


def test_top_n_larger_than_half_dedupes():
    p = _payload(top_n=5)
    ids = [row["feature_id"] for row in p["panels"][0]["rows"]]
    assert len(ids) == len(set(ids)) == 6


def test_all_features_are_aligned_by_the_first_pc():
    p = _payload(top_n=None)
    first = [row["feature_id"] for row in p["panels"][0]["rows"]]
    second = [row["feature_id"] for row in p["panels"][1]["rows"]]
    assert p["aligned"] is True and first == second
    values = [row["value"] for row in p["panels"][0]["rows"]]
    assert values == sorted(values, reverse=True)


def test_all_features_over_the_limit_is_rejected(monkeypatch):
    monkeypatch.setattr(pl, "MAX_ALL_FEATURES", 5)
    with pytest.raises(ValueError):
        _payload(top_n=None)


def test_r_requested_but_unavailable_falls_back_to_coefficient():
    p = _payload(r_available=False)
    assert p["value"] == "coefficient" and p["value_requested"] == "r"
    assert any("coefficient" in c for c in p["caveats"])


@pytest.mark.parametrize("pcs", [[3], [0], [1, 1], [1, 2, 3, 4]])
def test_invalid_pcs_are_rejected(pcs):
    with pytest.raises(ValueError):
        _payload(pcs=pcs)


def test_render_one_panel_per_pc_with_class_legend():
    fig = pl.render_loadings_plot(_payload(top_n=None))
    assert len([ax for ax in fig.axes if ax.get_visible()]) == 2
    plt.close(fig)
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_loadings_plot.py -q`
Expected: FAIL（`ImportError`）

- [ ] **Step 3: `metabolomix/plots/pca_loadings.py` を書く**

```python
"""PCA ローディング図（MCP 非依存）。ARF・分子種 PCA・mzTab の 3 経路が同じ形に変換して渡す。

1 パネル = 1 主成分の横棒。値は相関 r（既定。−1〜1、autoscale した PCA でだけ出せる）か
固有ベクトルの成分（coefficient）。top_n は主成分ごとに正の上位 N と負の上位 N（重複は 1 回）。
top_n=None は全件で、全パネルの行を最初の主成分の値の順にそろえる（行の並びで見比べられるように）。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §6。
"""
from __future__ import annotations

LOADINGS_SCHEMA = "lipidmix.pca_loadings.v1"
MAX_ALL_FEATURES = 60
MAX_PCS = 3
_CLASS_PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7e57c2",
                  "#00897b", "#8d6e63", "#5c6bc0", "#c0ca33"]


def _check_pcs(pcs, n_available):
    if not pcs or len(pcs) > MAX_PCS or len(set(pcs)) != len(pcs):
        raise ValueError(f"pcs は重複のない 1〜{MAX_PCS} 個の主成分番号で指定してください（受け取った値: {pcs!r}）。")
    for pc in pcs:
        if isinstance(pc, bool) or not isinstance(pc, int) or not 1 <= pc <= n_available:
            raise ValueError(f"pcs の {pc!r} はありません（この PCA の主成分は 1〜{n_available}）。")


def _row(feature, k, use_r):
    coef = float(feature["coefficient"][k])
    r = float(feature["r"][k]) if feature.get("r") is not None else None
    return {"feature_id": feature["feature_id"], "label": feature["label"], "ontology": feature.get("ontology"),
            "adduct": feature.get("adduct") or "", "mz": feature.get("mz"), "rt": feature.get("rt"),
            "coefficient": round(coef, 4), "r": round(r, 4) if r is not None else None,
            "value": round(r if use_r else coef, 4)}


def build_loadings_payload(features, *, source, result_id, explained_variance_ratio, pcs, top_n, value,
                           r_available, caveats=None):
    """特徴量の一覧から、主成分ごとに描く行を選んで payload にする。"""
    if value not in ("r", "coefficient"):
        raise ValueError(f"value は r / coefficient のどちらかを指定してください（受け取った値: {value!r}）。")
    if top_n is not None and (isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1):
        raise ValueError(f"top_n は 1 以上の整数か None で指定してください（受け取った値: {top_n!r}）。")
    n_available = min(len(explained_variance_ratio), min((len(f["coefficient"]) for f in features), default=0))
    _check_pcs(pcs, n_available)
    caveats = list(caveats or [])
    use_r = value == "r" and r_available
    if value == "r" and not r_available:
        caveats.append("この PCA は相関 r を出せない（autoscale していない、または特異値が無い）ため、"
                       "固有ベクトルの成分（coefficient）で描きました。")
    aligned = top_n is None
    if aligned and len(features) > MAX_ALL_FEATURES:
        raise ValueError(f"特徴量が {len(features)} あり、全件は描けません（{MAX_ALL_FEATURES} まで）。"
                         "top_n で絞ってください。")
    panels = []
    order = None
    for pc in pcs:
        k = pc - 1
        rows = [_row(f, k, use_r) for f in features]
        if aligned:
            if order is None:
                order = [row["feature_id"] for row in sorted(rows, key=lambda r: r["value"], reverse=True)]
            by_id = {row["feature_id"]: row for row in rows}
            chosen = [by_id[fid] for fid in order]
        else:
            ranked = sorted(rows, key=lambda r: r["value"], reverse=True)
            picked = ranked[:top_n] + ranked[-top_n:]
            seen = set()
            chosen = []
            for row in sorted(picked, key=lambda r: r["value"], reverse=True):
                if row["feature_id"] not in seen:
                    seen.add(row["feature_id"])
                    chosen.append(row)
        panels.append({"pc": f"PC{pc}", "explained_pct": round(float(explained_variance_ratio[k]) * 100, 1),
                       "rows": chosen})
    return {"plot_schema": LOADINGS_SCHEMA, "source": source, "result_id": result_id,
            "value": "r" if use_r else "coefficient", "value_requested": value, "top_n": top_n,
            "aligned": aligned, "panels": panels, "caveats": caveats}


def render_loadings_plot(payload, *, title=None):
    """payload を主成分ごとの横棒にする（色はクラス）。"""
    import matplotlib.pyplot as plt

    panels = payload["panels"]
    if not panels or not any(p["rows"] for p in panels):
        raise ValueError("描くローディングがありません。")
    classes = sorted({row["ontology"] for p in panels for row in p["rows"] if row.get("ontology")})
    color_of = {c: _CLASS_PALETTE[i % len(_CLASS_PALETTE)] for i, c in enumerate(classes)}
    n_rows = max(len(p["rows"]) for p in panels)
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 0.19 * n_rows + 1.6),
                             sharey=payload["aligned"], squeeze=False)
    is_r = payload["value"] == "r"
    for ax, panel in zip(axes[0], panels):
        rows = list(reversed(panel["rows"]))                 # 上ほど値が大きい
        ys = range(len(rows))
        ax.barh(list(ys), [row["value"] for row in rows], height=0.72,
                color=[color_of.get(row.get("ontology"), "#9e9e9e") for row in rows])
        ax.axvline(0, color="#6b6a64", lw=0.8)
        if is_r:
            ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel(f"{'loading r' if is_r else 'coefficient'} ({panel['pc']}, {panel['explained_pct']}%)",
                      fontsize=8)
        if not payload["aligned"] or ax is axes[0][0]:
            ax.set_yticks(list(ys))
            ax.set_yticklabels([f"{row['label']} {row['adduct']}".strip() for row in rows], fontsize=6.5)
        ax.grid(axis="x", color="#e4e3dc", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    if classes:
        handles = [plt.Rectangle((0, 0), 1, 1, color=color_of[c]) for c in classes]
        fig.legend(handles, classes, loc="lower center", ncol=min(len(classes), 8), frameon=False, fontsize=7)
    default = f"PCA loadings ({payload['source']}; {'r = correlation with the PC score' if is_r else 'eigenvector coefficient'})"
    fig.suptitle(title or default, fontsize=9)
    fig.tight_layout(rect=(0, 0.45 / fig.get_figheight() if classes else 0, 1, 0.97))
    return fig
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_loadings_plot.py -q`
Expected: PASS

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/plots/pca_loadings.py tests/test_pca_loadings_plot.py
git commit -m "feat(plots): PCA ローディング図の組み立てと描画" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t9.log 2>&1
```

---

### Task 10: `plot_pca_loadings` ツールと `save_figure(kind="pca_loadings")`

**Files:**
- Create: `metabolomix/tools/pca_loadings_tools.py`
- Modify: `server.py`、`metabolomix/core/session_state.py`（`AnalysisSession.last_loadings_plot`）、`metabolomix/tools/reports.py`
- Test: `tests/test_pca_loadings_tools.py`（新規）、`tests/test_save_figure.py`、`tests/test_missing_state_contract.py`
- 台帳 7 か所（`docs/workflow/plots.md`）、`docs/output_format/arf.md`、`docs/output_format/mztab.md`

**Interfaces:**
- Consumes: `last_pca_plot` の追加キー（Task 8）、`last_species_pca.loadings_rows`（Task 7）、`ds.last_pca`（`loadings` `singular_values` `explained_variance_ratio` `n_samples` `provenance`）と `ds.pp_feature_names` `ds.feature_annotations`、`get_pca_loading_features`（`arf/reader.py`）、`build_loadings_payload` / `render_loadings_plot`（Task 9）、`select_result`（`plots/result_output.py`）
- Produces: `pca_loadings_tools.plot_pca_loadings(source="auto", result_id=None, pcs=None, top_n=15, value="r", title=None, output=None)`、`session.last_loadings_plot = {"payload", "title"}`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_pca_loadings_tools.py`:

```python
"""plot_pca_loadings（ツール層）: 3 経路の PCA のローディング図。"""
from __future__ import annotations

import json
import math

import pytest

import server
from metabolomix.core import session_state


def _arf_plot(n_feat=4, with_loadings=True):
    plot = {"title": "PCA", "x_label": "PC1", "y_label": "PC2",
            "points": [{"x": 1.0, "y": 0.0, "label": "a"}, {"x": -1.0, "y": 0.5, "label": "b"},
                       {"x": 0.0, "y": -0.5, "label": "c"}, {"x": 0.3, "y": 0.2, "label": "d"}],
            "provenance": {"result_id": "arf-1", "dataset_id": "d"}}
    if with_loadings:
        plot.update({"loadings": [[0.5, -0.5, 0.5, -0.5][:n_feat], [0.1, 0.2, -0.3, 0.9][:n_feat]],
                     "singular_values": [4.0, 2.0], "explained_variance_ratio": [0.7, 0.2],
                     "feature_names": [f"Spot_{i}_height" for i in range(n_feat)], "n_fit": 4,
                     "scaling": "autoscale"})
    return plot


@pytest.fixture(autouse=True)
def _fresh():
    session_state.session = server.AnalysisSession()
    session_state.session.arf.features = [{"MasterAlignmentID": i, "Name": f"PG {30 + i}:0",
                                           "MassCenter": 700.0 + i, "RT": 3.0} for i in range(4)]
    yield


def _call(**kw):
    from metabolomix.tools.pca_loadings_tools import plot_pca_loadings
    return json.loads(plot_pca_loadings(output="payload", **kw))


def test_no_pca_is_missing_state():
    err = _call()["error"]
    assert err["code"] == "missing_state" and "arf_pca_species" in err["required_tools"]


def test_old_arf_result_without_loadings_is_not_a_candidate():
    session_state.session.arf.last_pca_plot = _arf_plot(with_loadings=False)
    assert _call()["error"]["code"] == "missing_state"


def test_arf_source_uses_reader_selection_and_computes_r():
    session_state.session.arf.last_pca_plot = _arf_plot()
    p = _call(top_n=1)
    assert p["source"] == "arf" and p["value"] == "r" and p["result_id"] == "arf-1"
    rows = {row["feature_id"]: row for row in p["panels"][0]["rows"]}
    assert set(rows) == {"0", "1"} or len(rows) == 2
    some = next(iter(rows.values()))
    assert some["r"] == pytest.approx(some["coefficient"] * 4.0 / math.sqrt(4), abs=1e-3)
    assert some["label"].startswith("PG")


def test_species_source_uses_stored_rows():
    session_state.session.arf.last_species_pca = {
        "points": [{"x": 0, "y": 0, "label": "a"}], "explained_variance_ratio": [0.8, 0.1],
        "loadings_rows": [{"feature_id": "7", "label": "DGDG 16:0_18:1", "ontology": "DGDG", "adduct": "[M+CH3COO]-",
                           "mz": 977.6, "rt": 4.1, "coefficient": [0.7, 0.1], "r": [0.99, 0.05]},
                          {"feature_id": "8", "label": "PG 16:0_19:1", "ontology": "PG", "adduct": "[M-H]-",
                           "mz": 761.5, "rt": 3.7, "coefficient": [-0.7, 0.1], "r": [-0.98, 0.04]}],
        "provenance": {"result_id": "sp-1", "dataset_id": "d"}}
    p = _call(top_n=None)
    assert p["source"] == "species" and p["aligned"] is True
    assert [row["feature_id"] for row in p["panels"][0]["rows"]] == ["7", "8"]


def test_two_sources_need_a_choice_and_bad_pcs_is_an_error():
    session_state.session.arf.last_pca_plot = _arf_plot()
    session_state.session.arf.last_species_pca = {
        "points": [{"x": 0, "y": 0, "label": "a"}], "explained_variance_ratio": [0.8, 0.1],
        "loadings_rows": [{"feature_id": "7", "label": "x", "ontology": None, "adduct": "", "mz": None, "rt": None,
                           "coefficient": [0.7, 0.1], "r": [0.9, 0.1]},
                          {"feature_id": "8", "label": "y", "ontology": None, "adduct": "", "mz": None, "rt": None,
                           "coefficient": [-0.7, 0.1], "r": [-0.9, 0.1]}],
        "provenance": {"result_id": "sp-1", "dataset_id": "d"}}
    assert _call()["error"]["code"] == "AMBIGUOUS_RESULT_SOURCE"
    out = _call(source="arf", pcs=[4])
    assert out["status"] == "error" and "pcs" in out["message"]


def test_image_and_save_figure(tmp_path, monkeypatch):
    from mcp.server.fastmcp import Image
    from metabolomix.tools.pca_loadings_tools import plot_pca_loadings
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    session_state.session.arf.last_pca_plot = _arf_plot()
    out = plot_pca_loadings()
    assert isinstance(out[1], Image) and "source=arf" in out[0]
    msg = server.save_figure("pca_loadings", "L1")
    png = tmp_path / "reports" / "figures" / "L1_pca_loadings.png"
    assert png.is_file() and png.with_suffix(".svg").is_file() and str(png) in msg


def test_failed_call_clears_previous_figure():
    session_state.session.arf.last_pca_plot = _arf_plot()
    _call()
    assert session_state.session.last_loadings_plot is not None
    _call(value="loading")
    assert session_state.session.last_loadings_plot is None
```

mzTab 経路のテストも同じファイルに足す（`tests/test_dataset_analysis_tools.py` の `_load_ds` を使う）:

```python
def test_mztab_source_draws_dataset_pca():
    from tests.test_dataset_analysis_tools import _load_ds
    from metabolomix.tools.dataset_analysis_tools import dataset_pca, dataset_preprocess
    _load_ds()
    dataset_preprocess()
    dataset_pca(n_components=2)
    p = _call(source="mztab", top_n=3)
    assert p["source"] == "mztab" and p["value"] == "r"
    assert p["panels"][0]["rows"][0]["label"].startswith("Compound")
```

`tests/test_save_figure.py` の parametrize に `("pca_loadings", "pca_loadings_plot", ["plot_pca_loadings"])` と `"pca_loadings"`（source を拒否する側）を足す。`tests/test_missing_state_contract.py` に:

```python
    def test_plot_pca_loadings(self):
        self.assert_missing(
            server.plot_pca_loadings(output="payload"),
            "pca_result",
            ["arf_parser", "arf_pca_preprocessed", "arf_pca_species", "dataset_pca"],
            "PCA",
        )
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_loadings_tools.py -q`
Expected: FAIL（`ModuleNotFoundError`）

- [ ] **Step 3: セッションの欄を足す**

`metabolomix/core/session_state.py` の `AnalysisSession.__init__` に `self.last_loadings_plot = None  # plot_pca_loadings / save_figure(kind="pca_loadings") が参照（どの経路の PCA でも 1 か所）` を足す。

- [ ] **Step 4: `metabolomix/tools/pca_loadings_tools.py` を書く**

```python
"""PCA のローディング図（plot_pca_loadings）。ARF・分子種 PCA・mzTab の 3 経路を同じ図にする。

どの PCA を描くかは save_figure(kind="pca") と同じ規則で選ぶ（どれも優先せず、曖昧なら止める）。
ARF は上位 N の選び方とスポット情報に既存の arf/reader.py get_pca_loading_features() を使う。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §6。
"""
from __future__ import annotations

import math

from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.core import mcp_errors, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.plots import pca_loadings as loadings_plot
from metabolomix.plots import render as plot_render

__all__ = ["plot_pca_loadings"]


def _candidates() -> list[dict]:
    from metabolomix.analysis.result_state import is_current

    out = []
    arf = session_state.session.arf
    plot = getattr(arf, "last_pca_plot", None)
    if plot and plot.get("loadings") and plot.get("feature_names"):
        prov = plot.get("provenance") or {}
        out.append({"source": "arf", "result_id": prov.get("result_id"), "dataset_id": prov.get("dataset_id"),
                    "valid": True, "result": plot})
    species = getattr(arf, "last_species_pca", None)
    if species and species.get("loadings_rows"):
        prov = species.get("provenance") or {}
        out.append({"source": "species", "result_id": prov.get("result_id"), "dataset_id": prov.get("dataset_id"),
                    "valid": True, "result": species})
    ds = getattr(session_state.session, "dataset", None)
    ds_pca = getattr(ds, "last_pca", None) if ds is not None else None
    if ds_pca and ds_pca.get("loadings"):
        prov = ds_pca.get("provenance") or {}
        out.append({"source": "mztab", "result_id": prov.get("result_id"), "dataset_id": prov.get("dataset_id"),
                    "valid": is_current(ds, ds_pca), "result": ds_pca, "ds": ds})
    return out


def _arf_features(plot: dict, pcs: list[int], top_n: int | None) -> tuple[list[dict], bool]:
    from metabolomix.arf.reader import get_pca_loading_features
    from metabolomix.arf.selection import display_name
    from metabolomix.arf.tools import _sibling_arf2_path
    from metabolomix.arf2 import reader as arf2_reader

    names = plot["feature_names"]
    loadings = plot["loadings"]
    if max(pcs) > len(loadings):
        raise ValueError(f"pcs の {max(pcs)} はありません（この PCA の主成分は 1〜{len(loadings)}）。")
    if top_n is None and len(names) > loadings_plot.MAX_ALL_FEATURES:
        raise ValueError(f"特徴量が {len(names)} あり、全件は描けません（{loadings_plot.MAX_ALL_FEATURES} まで）。"
                         "top_n で絞ってください。")
    picked = get_pca_loading_features(
        {"loadings": loadings, "explained_variance_ratio": plot.get("explained_variance_ratio") or []},
        session_state.session.arf.features or [], names, top_n=top_n or len(names), n_pcs=max(pcs))
    arf2_path = _sibling_arf2_path()
    catalog = ({int(r["MasterAlignmentID"]): r for r in arf2_reader.load_catalog(str(arf2_path))}
               if arf2_path else {})
    column_of = {}
    for j, name in enumerate(names):
        parts = name.split("_")
        if len(parts) >= 3 and parts[1].isdigit():
            column_of.setdefault(int(parts[1]), j)
    ids = []
    for pc in picked:
        for item in pc["positive"] + pc["negative"]:
            if item["id"] not in ids:
                ids.append(item["id"])
    r_ok = plot.get("scaling") == "autoscale" and plot.get("singular_values") and plot.get("n_fit")
    factors = ([float(s) / math.sqrt(plot["n_fit"]) for s in plot["singular_values"]] if r_ok else None)
    features = []
    meta = {item["id"]: item for pc in picked for item in pc["positive"] + pc["negative"]}
    for sid in ids:
        j = column_of[sid]
        coef = [float(loadings[k][j]) for k in range(len(loadings))]
        row = catalog.get(sid, {})
        features.append({
            "feature_id": str(sid),
            "label": display_name(row.get("Name") or "") or meta[sid].get("annotation") or f"Spot {sid}",
            "ontology": (row.get("Ontology") or "").strip() or None, "adduct": row.get("AdductType") or "",
            "mz": meta[sid].get("m_z"), "rt": meta[sid].get("rt"), "coefficient": coef,
            "r": [c * factors[k] for k, c in enumerate(coef)] if factors else None})
    return features, bool(factors)


def _mztab_features(chosen: dict) -> tuple[list[dict], bool]:
    last = chosen["result"]
    ds = chosen["ds"]
    loadings = last["loadings"]
    s = last.get("singular_values")
    n = last.get("n_samples")
    factors = [float(v) / math.sqrt(n) for v in s] if s and n else None
    annotations = getattr(ds, "feature_annotations", {}) or {}
    metadata = getattr(ds, "feature_metadata", {}) or {}
    features = []
    for j, fid in enumerate(ds.pp_feature_names):
        coef = [float(loadings[k][j]) for k in range(len(loadings))]
        ann = annotations.get(fid) or {}
        meta = metadata.get(fid) or {}
        features.append({"feature_id": str(fid), "label": ann.get("name") or meta.get("name") or str(fid),
                         "ontology": None, "adduct": "", "mz": meta.get("mz"), "rt": meta.get("rt"),
                         "coefficient": coef,
                         "r": [c * factors[k] for k, c in enumerate(coef)] if factors else None})
    return features, bool(factors)


def _caption(payload: dict) -> str:
    parts = "、".join(f"{p['pc']} {p['explained_pct']}%（{len(p['rows'])} 行）" for p in payload["panels"])
    caption = (f"PCA ローディング（source={payload['source']} result_id={payload['result_id'] or '(なし)'}、"
               f"値={payload['value']}）: {parts}。")
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return caption


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def plot_pca_loadings(
    source: str = "auto",
    result_id: str | None = None,
    pcs: list[int] | None = None,
    top_n: int | None = 15,
    value: str = "r",
    title: str | None = None,
    output: str | None = None,
) -> list | str:
    """保存済みの PCA（arf_parser / arf_pca_preprocessed / arf_pca_species / dataset_pca）のローディングを横棒で描く。

    - source: "auto"（既定）/ "arf" / "species" / "mztab"。auto はどれも優先せず、候補が 2 つ以上あると
      AMBIGUOUS_RESULT_SOURCE で止まる。result_id で特定の結果を名指しできる。
    - pcs: 描く主成分（既定 [1, 2]、1〜3 個）。
    - top_n: 主成分ごとに正の上位 N と負の上位 N（既定 15）。None は全件（60 特徴量まで）で、
      全パネルの行を最初の主成分の順にそろえる。
    - value: "r"（既定。各特徴量と主成分スコアの相関、−1〜1）/ "coefficient"（固有ベクトルの成分）。
      r は autoscale した PCA でだけ出せ、出せないときは coefficient で描いて説明文に書く。
    - 棒の色はクラス（ARF・分子種 PCA は .arf2 の Ontology）。
    - output: "image"（既定）/ "payload"（lipidmix.pca_loadings.v1）。ファイルは save_figure(kind="pca_loadings")。
    """
    session = session_state.session
    session.last_loadings_plot = None
    try:
        mode = plot_render.resolve_plot_output(output)
        pcs = list(pcs) if pcs is not None else [1, 2]
        if value not in ("r", "coefficient"):
            raise ValueError(f"value は r / coefficient のどちらかを指定してください（受け取った値: {value!r}）。")
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})
    candidates = _candidates()
    if not candidates:
        return mcp_errors.missing_state(
            "pca_result", ["arf_parser", "arf_pca_preprocessed", "arf_pca_species", "dataset_pca"],
            "先に PCA を実行してください（ローディングを持つ PCA 結果がありません）。")
    from metabolomix.plots.result_output import select_result
    try:
        chosen = select_result(candidates, source=source, result_id=result_id)
    except DomainError as exc:
        return mcp_errors.mztab_error(exc.code, exc.message, exc.details or None)
    try:
        result = chosen["result"]
        if chosen["source"] == "arf":
            features, r_ok = _arf_features(result, pcs, top_n)
        elif chosen["source"] == "species":
            features, r_ok = list(result["loadings_rows"]), True
        else:
            features, r_ok = _mztab_features(chosen)
        payload = loadings_plot.build_loadings_payload(
            features, source=chosen["source"], result_id=chosen.get("result_id"),
            explained_variance_ratio=result.get("explained_variance_ratio") or [], pcs=pcs, top_n=top_n,
            value=value, r_available=r_ok)
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})
    if mode == plot_render.PAYLOAD:
        session.last_loadings_plot = {"payload": payload, "title": title}
        return json_payload(payload)
    png = plot_render.figure_to_png(loadings_plot.render_loadings_plot(payload, title=title))
    session.last_loadings_plot = {"payload": payload, "title": title}
    return [_caption(payload), Image(data=png, format="png")]
```

`server.py` の `from metabolomix.tools.dataset_analysis_tools import *` の次に `from metabolomix.tools.pca_loadings_tools import *  # plot_pca_loadings` を足す（mzTab の DatasetState を読むので dataset 系の後）。

- [ ] **Step 5: `save_figure` に `pca_loadings` を足す**

`metabolomix/tools/reports.py`: `FIGURE_KINDS = ("pca", "volcano", "eic", "group_intensity", "species", "pca_loadings")`。docstring に `"pca_loadings"（plot_pca_loadings。dpi 300 の PNG と同名 .svg）` を足し、分岐の末尾を:

```python
    if kind == "species":
        return _save_species(analysis_id, title)
    return _save_pca_loadings(analysis_id, title)


def _save_pca_loadings(analysis_id: str, title: str | None) -> str:
    from metabolomix.plots.pca_loadings import render_loadings_plot

    last = getattr(session_state.session, "last_loadings_plot", None)
    if not last or not last.get("payload"):
        return mcp_errors.missing_state(
            "pca_loadings_plot", ["plot_pca_loadings"],
            "先に plot_pca_loadings を実行してください（ローディング図がありません）。")
    out_path = _figures_dir() / f"{knowledge_store.make_slug(analysis_id)}_pca_loadings.png"
    _write_png_and_svg(render_loadings_plot(last["payload"], title=title or last.get("title")), out_path)
    payload = last["payload"]
    return (f"ローディング図を保存: {out_path}（同名の .svg も保存。source={payload['source']} "
            f"result_id={payload['result_id'] or '(なし)'}）\n"
            f"本文に ![PCA loadings](figures/{out_path.name}) で埋め込めます。")
```

- [ ] **Step 6: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pca_loadings_tools.py tests/test_save_figure.py tests/test_missing_state_contract.py -q`
Expected: PASS

- [ ] **Step 7: 台帳 7 か所を直す（73 → 74）**

1. `EXPECTED_TOOLS` に `"plot_pca_loadings"`、コメント `# PCA のローディング図: plot_pca_loadings を追加(73→74)。`
2. `"plot_pca_loadings": READ_ONLY,`
3. `"plots.md": ("save_figure", "plot_pca_loadings"),`、合計 `53`
4. `docs/workflow/plots.md` の payload 表に `| lipidmix.pca_loadings.v1 | plot_pca_loadings(output="payload") | save_figure(kind="pca_loadings") | 主成分ごとの横棒。値は r か coefficient。PNG（dpi 300）と SVG |` を足し、`## save_figure` の連鎖の末尾を:

```
12.├─ [kind=species] metabolomix/tools/reports.py  _save_species()
13.│  └─ metabolomix/plots/species.py  render_species_plot()
14.└─ [kind=pca_loadings] metabolomix/tools/reports.py  _save_pca_loadings()
15.   └─ metabolomix/plots/pca_loadings.py  render_loadings_plot()
```

   末尾に節を足す:

```markdown
## plot_pca_loadings

前提: ローディングを持つ PCA 結果がある（`arf_parser` / `arf_pca_preprocessed` / `arf_pca_species` /
`dataset_pca`。無ければ `MissingState`）
状態変更: `session.last_loadings_plot`（payload・title）を更新。ファイルは書かない。

入力元は 3 つ（手順 2）。選び方は `save_figure(kind="pca")` と同じで、どれも優先せず曖昧なら
`AMBIGUOUS_RESULT_SOURCE` で止まる（手順 3）。ARF は上位 N の選び方とスポット情報に
`get_pca_loading_features()` を使い（手順 4〜5）、分子種 PCA は保存済みの行を、mzTab は
`pp_feature_names` と特徴量の注釈を使う（手順 6）。r は autoscale した PCA でだけ
成分 × 特異値 / √n で出し、出せなければ coefficient で描く（手順 7）。

1. metabolomix/tools/pca_loadings_tools.py  plot_pca_loadings()
2. └─ metabolomix/tools/pca_loadings_tools.py  _candidates()
3. └─ metabolomix/plots/result_output.py  select_result()
4. ├─ [source=arf] metabolomix/tools/pca_loadings_tools.py  _arf_features()
5. │  └─ metabolomix/arf/reader.py  get_pca_loading_features()
6. ├─ [source=mztab] metabolomix/tools/pca_loadings_tools.py  _mztab_features()
7. └─ metabolomix/plots/pca_loadings.py  build_loadings_payload()
8. ├─ [output=payload] metabolomix/core/serialization.py  json_payload()
9. └─ [output=image] metabolomix/plots/pca_loadings.py  render_loadings_plot()
10.   └─ metabolomix/plots/render.py  figure_to_png()
```

5. `docs/workflow/index.md`: plots.md の行を `| [plots.md](plots.md) | 図の保存（save_figure）・PCA ローディング図と payload 契約の比較 | 2 |`、`合計 53 ツール`、`登録ツール総数は 74（53 + 21）`
6. `USAGE.md`: 見出し `(全74ツール)`、`save_figure` の行の前に:

   `| `plot_pca_loadings` | 保存済みの PCA(`arf_parser` / `arf_pca_preprocessed` / `arf_pca_species` / `dataset_pca`)のローディングを主成分ごとの横棒で描く(`source`, `result_id`, `pcs`, `top_n`, `value`, `title`, `output`)。`source` は `auto`/`arf`/`species`/`mztab` で、`auto` はどれも優先せず候補が 2 つ以上あると `AMBIGUOUS_RESULT_SOURCE`。`pcs` 既定 [1, 2](1〜3 個)。`top_n`(既定 15)は主成分ごとに正の上位 N と負の上位 N、`None` は全件(60 特徴量まで)で行を最初の主成分の順にそろえる。`value="r"`(既定。特徴量と主成分スコアの相関)は autoscale した PCA でだけ出せ、出せなければ `coefficient` で描いて説明文に書く。棒の色はクラス。ローディングを持たない古い ARF の PCA 結果は候補にならない(PCA をやり直す)。既定は PNG 画像と説明、`output="payload"` で `lipidmix.pca_loadings.v1`。ファイルは `save_figure(kind="pca_loadings")`。 |`

   `save_figure` の行の kind の列挙に `/ `pca_loadings`(plot_pca_loadings。dpi 300 の PNG と同名 `.svg`)` を足す。
7. `CLAUDE.md`: `ツール 74・リソース 5・リソーステンプレート 3`

`docs/output_format/arf.md` に `### 11.8 `plot_pca_loadings` の返り値（`lipidmix.pca_loadings.v1`）` を足す（`value` は実際に描いた値の種類で `value_requested` と違えば caveats に理由、`r` は −1〜1 で主成分の符号は任意、`aligned` の意味、ARF のローディングは前処理後の全スポットのうち上位だけ）。`docs/output_format/mztab.md` に `plot_pca_loadings(source="mztab")` が `dataset_pca` の結果を描けることを 1 段落で足す。

- [ ] **Step 8: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS

- [ ] **Step 9: コミット（バックグラウンド）**

```bash
git add -A metabolomix server.py tests docs/workflow docs/output_format USAGE.md CLAUDE.md
git commit -m "feat(pca): 3 経路の PCA のローディング図 plot_pca_loadings と save_figure(kind=pca_loadings) を足す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t10.log 2>&1
```

---

### Task 11: 実データでの確認と記録（20260930_EV）

**Files:**
- Create（リポジトリ外・scratchpad）: `verify_species_pca.py`
- Modify: `docs/HISTRY.md`（末尾に追記）、`docs/task.md`、`C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`、`metabolomix/core/mcp_core.py` の `MCP_INSTRUCTIONS`（保存ツールの案内があれば）

- [ ] **Step 1: 試作（v3）と同じ数値になるかを確かめるスクリプトを scratchpad に書く**

```python
"""実データ（20260930_EV NEG）で新ツールが v3 の試作と同じ数値を出すかを確かめる（読み取りのみ）。"""
import json
import sys

sys.path.insert(0, r"C:\Users\yuu18\Metabolomix_with_LLM")
import server  # noqa: F401
from metabolomix.arf.species_tools import arf_plot_species, arf_pca_species
from metabolomix.arf.tools import arf_exclude, arf_parser
from metabolomix.core import session_state
from metabolomix.tools.pca_loadings_tools import plot_pca_loadings

arf_parser(file_path=r"C:\Users\yuu18\datasets\20260930_EV\neg_wiff\AlignResult-20261051221_PeakProperties.arf",
           top_features=1)
groups = ["EV_cell", "EV_PlnA", "cell_C", "cell_PlnA"]
common = dict(groups=groups, low_reliability_samples=["20260930_cell_PlnA_NEG_1"],
              exclude_auto_likely_wrong=True, standard_samples=["systemlot8"], require_msms=True)
share = json.loads(arf_plot_species(items=["DGDG"], share_basis=["PG", "DGDG", "MGDG", "CL"],
                                    output="payload", **common))
row = next(s for s in share["spots"] if s["label"] == "DGDG 18:1_18:1")
means = {g["label"]: g["mean"] for g in row["groups"]}
print("DGDG 18:1_18:1 share", means)        # v3: cell_C 1.43, cell_PlnA 2.42, EV_cell 2.49
pca = json.loads(arf_pca_species(items=["PG", "DGDG", "MGDG", "CL"], groups=["cell_C", "cell_PlnA"],
                                 low_reliability_samples=["20260930_cell_PlnA_NEG_1"],
                                 exclude_auto_likely_wrong=True, standard_samples=["systemlot8"],
                                 require_msms=True, normalize="total", orient_by="cell_PlnA", output="payload"))
print("PCA evr", pca["explained_variance_ratio"][:2])   # v3: PC1 82.7 %, PC2 11.1 %
print("PC1 r top", pca["pc1_r_top"][:3])                 # v3: DGDG 18:1_18:1 +1.00, 18:1_19:1 +1.00, 18:0_18:1 +1.00
ld = json.loads(plot_pca_loadings(source="species", top_n=None, output="payload"))
print("loadings rows", len(ld["panels"][0]["rows"]))     # v3: 38
```

Run: `C:/Python314/python.exe -I <scratchpad>/verify_species_pca.py`
Expected: 割合の平均がスライドの値（cell_C 1.43、cell_PlnA 2.42、EV_cell 2.49）と小数第 2 位まで一致、PC1 82.7 %・PC2 11.1 %、PC1 の r の上位が 18:1 を持つ DGDG、ローディングの行数 38。ずれたら原因（採用の規則・群の解決・正規化）を調べて報告する（数値を合わせるために実装を曲げない）。

- [ ] **Step 2: 記録を書く**

- `docs/HISTRY.md` の末尾に `## 2026-10-08 分子種ごとの図・分子種 PCA・ローディング図・save_figure` の節を足す（経緯: 上司の要望 → v3 試作 → MCP 化、決まったこと: B 案で旧保存ツールを同時撤去・項目解決の共通化・r の定義、実データ確認の結果、Use-LLLM 側の変更）。
- `docs/task.md` は並列エージェントがいなければ、この作業で閉じた項目を消し、残課題（mzTab 経路の分子種ごとの図は対象外、など必要なもの）を足す。
- vault の流れ図に、ARF 経路の「分子種ごとの図」「分子種 PCA → ローディング図」の枝と、保存の一本化（`save_figure`）を書き、末尾の出典行の日付を 2026-10-08 にする。
- `grep -n "save_" metabolomix/core/mcp_core.py` で `MCP_INSTRUCTIONS` に旧ツール名や `save_*_figure` が無いことを確かめる（あれば `save_figure` に直し、全テストを回す）。

- [ ] **Step 3: コミット（リポジトリ内の変更があれば。バックグラウンド）**

```bash
git add metabolomix/core/mcp_core.py
git commit -m "docs: MCP_INSTRUCTIONS の図の保存の案内を save_figure に直す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > <scratchpad>/commit_t11.log 2>&1
```

（`docs/HISTRY.md` `docs/task.md` は追跡外なのでコミットしない。`MCP_INSTRUCTIONS` に変更が無ければこのコミットは作らない。）

---

### Task 12: Use-LLLM 側の追随（別リポジトリ）

**Files（`C:\Users\yuu18\Use-LLLM`）:**
- Modify: `tests/test_tool_catalog_completeness.py`、`tests/test_approval_snapshot.py`、`tests/test_general_agent.py`、`README.md`

- [ ] **Step 1: 作業ブランチを切る**

Use-LLLM には CLAUDE.md が無い。テストは README のとおり `uv run pytest -q` と `uv run ruff check src tests`。

```bash
cd C:/Users/yuu18/Use-LLLM && git status --short && git checkout -b chore/ms-data-parser-save-figure
```

- [ ] **Step 2: ツール名を直す**

- `tests/test_tool_catalog_completeness.py` の `MS_DATA_PARSER_TOOLS`: `"save_pca_figure"` `"save_volcano_figure"` `"save_eic_figure"` を消し、`"save_figure"` を足す（新しい描画ツール `"arf_plot_species"` `"arf_pca_species"` `"plot_pca_loadings"` もサーバのカタログに載るので足す）。
- `tests/test_approval_snapshot.py`: 3 行を消して `"save_figure": LOCAL_WRITE,` を足す。新しい描画ツールは `READ_ONLY` で足す。
- `tests/test_general_agent.py`: `"lipidmix::save_pca_figure"` を `"lipidmix::save_figure"` に（3 か所）。
- `README.md` 22 行目: `save_pca_figure`、`save_volcano_figure`、`save_eic_figure` の列挙を `save_figure` に直す。

- [ ] **Step 3: Use-LLLM のテストを回す**

Run: `cd C:/Users/yuu18/Use-LLLM && uv run pytest -q && uv run ruff check src tests`
Expected: 全件 PASS。ツール名の一覧やスナップショットが他のテストにもあって落ちたら、同じように直す。

- [ ] **Step 4: コミット（push はしない。ユーザーに確認する）**

```bash
git add tests/test_tool_catalog_completeness.py tests/test_approval_snapshot.py tests/test_general_agent.py README.md
git commit -m "chore: ms-data-parser の図の保存ツールの一本化（save_figure）に追随する" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## 完了後

- 本体リポジトリで `superpowers:finishing-a-development-branch` に従う。`main` へのマージは `--no-ff` で `Merge feat/species-plot-and-pca-loadings: 分子種ごとの図・分子種 PCA・ローディング図を足し、図の保存を save_figure に一本化する`。push は都度確認する。
- Use-LLLM のブランチのマージと push もユーザーに確認する。
