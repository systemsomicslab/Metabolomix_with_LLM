# MSLipidMapper パスウェイ投影ツール Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** ARF / DatasetState の前処理済み脂質データを MSLipidMapper（ローカル Rscript）でパスウェイへ投影し、LLM 向けのノード表・逆向き隣接ペアと、チャット表示用の PNG を返す MCP ツール `pathway_project` / `pathway_plot` を作る。

**Architecture:** Python 側が上流 `write_mztab_m` と同形の最小 mzTab-M を書き、`metabolomix/pathway/mslm_bridge.R` を `Rscript` で起動して上流の公開関数（`mztab_pathway_measurements` / `render_pathway_from_config`）を呼ぶ。返ってきた JSON-LD を Python で読み、既存の Welch＋BH で 2 群比較し、接合ノード経由の隣接で向きを判定する。PDF は PyMuPDF で PNG にして `Image` で返す。

**Tech Stack:** Python 3.14（`C:/Python314/python.exe`）、FastMCP、numpy、PyMuPDF（`fitz`）、R 4.6（上流 MSLipidMapper をローカルクローンからパッケージとしてインストール）。

**Spec:** `docs/superpowers/specs/2026-10-14-mslipidmapper-pathway-design.md`

## Global Constraints

- テストは `C:/Python314/python.exe -m pytest tests -q`（リポジトリルート）。R が無くても全部緑。
- 全ツールに `structured_output=False` と `ToolAnnotations`。JSON は `metabolomix.core.serialization.json_payload()`、float は丸める。
- 外部契約は旧名: 新しい環境変数は既存規則に合わせる（上流と共有する `MSLIPIDMAPPER_*` はそのまま上流の名前を使う）。設定ファイルは `lipidmix.local.toml` の `[mslipidmapper]`。読むのは `core/user_config.get_setting` だけ。
- `metabolomix.core.mcp_core` は leaf。`metabolomix/pathway/` の純ロジック（handoff / jsonld / compare / bridge）は `session_state` も `mcp_core` も import しない。import してよいのは `pathway/tools.py` だけ。
- 成果物は `mcp_core._resolve_report_dir() / "pathway" / <run_id>`。リポジトリ側に書かない。
- 前提状態が無いときは `mcp_errors.missing_state(...)`。確定エラーは新設 `mcp_errors.pathway_error(code, message, details)`。
- 差次的エクスポート（`analysis/export_contract.py`）は触らない。`CONTRACT_VERSION` 据え置き。
- 研究室参照ライブラリのパスや上流クローンの絶対パスを追跡文書・テストに書かない。
- fixture はテスト自身が作る（上流の例ファイルや `data/` に依存しない）。
- コミットはフックで全テストが走る（数分）。`git commit` はバックグラウンド実行し、出力をファイルへ落とす。実行中に作業ツリーを編集しない。`git stash` 禁止。
- コミットメッセージ末尾: `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>`

## Review Focus

1. 同じ ARF に `height` と `area` の 2 種の prop が載った前処理行列 → 1 スポットが 2 行として二重計上されず、レシピ先頭の prop だけが使われる（Task 2 にテスト）。
2. 群指定が因子トークンでプールされる（`group_a="PlnA"` が複数 Class ID を掴む）→ mzTab 上の群ラベルは比較の spec 文字列そのもので、PDF の x 軸も比較 2 群が並ぶ（Task 2 にテスト）。
3. JSON-LD の接合ノードが 3 本以上の枝を持つ（分岐）→ 枝の先の脂質ノード全ての組が隣接として扱われ、`PathwayEntity` を挟む先は含まれない（Task 4 にテスト）。
4. 片群の全サンプルでスコープ値が `null`（MissingValue）→ そのノードは `call=NA`、BH の分母に入らず、逆向き判定では判定不能に数える（Task 5 にテスト）。
5. Rscript がタイムアウトしたとき子孫（ヘッドレス Edge）まで残らない → `terminate_tree` が呼ばれ `MSLM_TIMEOUT` が返る（Task 6 にテスト）。

---

### Task 0: spec の訂正（実装前に確定した事実の反映）

**Files:**
- Modify: `docs/superpowers/specs/2026-10-14-mslipidmapper-pathway-design.md`

計画作成中のコード読解で判明した 3 点を spec に反映する。

- [ ] **Step 1: §2「配置」の上流の扱いを書き換える**

「上流パッケージはインストールしない。クローン内の R ファイルを `source()` する…」の段落を次で置き換える:

```markdown
- 上流はローカルクローンから **R パッケージとしてインストールする**
  （`R CMD INSTALL <clone>`。上流 REST ワーカーも `library(MSLipidMapper)` で動く）。
  上流の R ファイルの多くがトップレベルで `library()` を呼ぶため、PDF 描画に要る関数群
  （`aggregate_to_class_se` `get_plot_fun` `.render_pathway_cytoscape_pdf` ほか）を
  部分的に `source()` することはできない。ブリッジは `MSLipidMapper::` の公開関数だけを呼ぶ。
  設定は `[mslipidmapper] rscript`（必須）・`browser`（描画時）・`repo`（任意。来歴の
  コミット記録とインストール案内に使う）。環境変数 `MSLIPIDMAPPER_RSCRIPT`
  `MSLIPIDMAPPER_BROWSER` `MSLIPIDMAPPER_REPO`。ブラウザは描画 YAML の `pathway.browser` で渡す。
```

- [ ] **Step 2: §4「特徴（行）」のサブクラスの出どころを直す**

「arf はスポットの `Ontology`」を次で置き換える:

```markdown
- サブクラス: arf は**同じアラインメントの兄弟 `.arf2`** のスポット（`MasterAlignmentID` で結合）の
  `Ontology` と `Name`（ARF のスポットは Ontology を持たない）。兄弟 `.arf2` が無ければ
  `missing_state("sibling_arf2", ["load_dataset"])`（既存 `_require_arf_with_arf2` と同じ）。
  前処理行列が複数の prop（`height` と `area` 等）を持つときはレシピ `props` の先頭だけを使う。
```

- [ ] **Step 3: §4「サンプル（列）」の群の規則を具体化する**

「群: arf は…、dataset は…」を次で置き換える:

```markdown
- 群: 比較に入るサンプルは**比較の指定名**（arf は `group_a`/`group_b` の文字列、dataset は
  `group_a_label`/`group_b_label`）を群ラベルにする（因子トークンでプールした群が図でも 1 群に並ぶ）。
  比較に入らない生体試料は、arf は `sample_meta[...]["group"]`、dataset は実験情報シート
  （`sample_metadata_rows` の `group`、`include` が偽のものは除く）があればその群で書き、
  無ければ書かない。
```

- [ ] **Step 4: §3 の `pathway_project` 引数に dataset 用ラベルを足す**

シグネチャのコードブロックの `group_a=..., group_b=...,` の次の行に追加:

```
  group_a_label="group_a", group_b_label="group_b",   # dataset のときだけ使う（arf は指定文字列）
```

- [ ] **Step 5: コミット**（バックグラウンド、フック完走を待つ）

```bash
git add docs/superpowers/specs/2026-10-14-mslipidmapper-pathway-design.md
git commit -m "docs(spec): パスウェイ投影の spec を実装前の調査結果で訂正（上流のインストール・ARF の Ontology 出どころ・群規則）" > <scratchpad>/commit.log 2>&1
```

---

### Task 1: 設定キーとエラー封筒

**Files:**
- Modify: `metabolomix/core/user_config.py`（`SETTINGS`）
- Modify: `lipidmix.example.toml`
- Modify: `metabolomix/core/mcp_errors.py`
- Test: `tests/test_pathway_config.py`（新規）

**Interfaces:**
- Produces: 設定キー `"mslipidmapper.rscript"` `"mslipidmapper.browser"` `"mslipidmapper.repo"`、`mcp_errors.PATHWAY_ERROR_CODES: frozenset[str]`、`mcp_errors.pathway_error(code: str, message: str, details: dict | None = None) -> str`。

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pathway_config.py
import json

import pytest

from metabolomix.core import mcp_errors, user_config


def test_mslipidmapper_settings_read_env_first(monkeypatch, tmp_path):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[mslipidmapper]\nrscript = 'R/bin/Rscript.exe'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    got = user_config.get_setting("mslipidmapper.rscript")
    assert got.source == "config_file"
    assert got.value == str(tmp_path / "R" / "bin" / "Rscript.exe")
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", "X:/Rscript.exe")
    assert user_config.get_setting("mslipidmapper.rscript").value == "X:/Rscript.exe"


@pytest.mark.parametrize("key,env", [
    ("mslipidmapper.rscript", "MSLIPIDMAPPER_RSCRIPT"),
    ("mslipidmapper.browser", "MSLIPIDMAPPER_BROWSER"),
    ("mslipidmapper.repo", "MSLIPIDMAPPER_REPO"),
])
def test_mslipidmapper_setting_env_names(key, env):
    assert user_config.SETTINGS[key] == env


def test_pathway_error_envelope():
    body = json.loads(mcp_errors.pathway_error("MSLM_TIMEOUT", "遅い", {"timeout_s": 300}))
    assert body == {"error": {"code": "MSLM_TIMEOUT", "message": "遅い",
                              "details": {"timeout_s": 300}}}


def test_pathway_error_rejects_unknown_code():
    with pytest.raises(ValueError):
        mcp_errors.pathway_error("NOPE", "x")


def test_pathway_error_codes_are_the_spec_set():
    assert mcp_errors.PATHWAY_ERROR_CODES == frozenset({
        "MSLM_NOT_CONFIGURED", "MSLM_R_PACKAGES_MISSING", "MSLM_RUN_FAILED",
        "MSLM_TIMEOUT", "MSLM_OUTPUT_INVALID", "PATHWAY_RENDER_FAILED",
        "PATHWAY_BAD_REQUEST"})
```

`tests/conftest.py` は `user_config.SETTINGS` の全環境変数を自動で消し、`LIPIDMIX_CONFIG` を存在しないパスへ向けている。`SETTINGS` に足せば新しい 3 変数も自動で消えるので conftest は触らない。

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_config.py -q`
Expected: FAIL（`KeyError: 'mslipidmapper.rscript'` / `AttributeError: pathway_error`）

- [ ] **Step 3: 実装**

`metabolomix/core/user_config.py` の `SETTINGS` に追加:

```python
    "mslipidmapper.rscript": "MSLIPIDMAPPER_RSCRIPT",
    "mslipidmapper.browser": "MSLIPIDMAPPER_BROWSER",
    "mslipidmapper.repo": "MSLIPIDMAPPER_REPO",
```

`lipidmix.example.toml` の末尾に追加:

```toml

[mslipidmapper]
# pathway_project / pathway_plot が呼ぶ Rscript。MSLipidMapper を R パッケージとして
# インストール済みであること（R CMD INSTALL <クローン>）。
# rscript = 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe'
# pathway_plot の PDF 描画に使うヘッドレス Chromium（Chrome / Edge）。
# browser = 'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
# 上流クローンの場所。任意（来歴に上流のコミットを記録し、インストール案内に使う）。
# repo = 'C:\path\to\MSLipidMapper'
```

`metabolomix/core/mcp_errors.py` の `mztab_error` の後に追加:

```python
PATHWAY_ERROR_CODES = frozenset({
    # [mslipidmapper] rscript が未設定、または指す先が無い。
    "MSLM_NOT_CONFIGURED",
    # MSLipidMapper（R パッケージ）か、その公開関数が無い。自動インストールはしない。
    "MSLM_R_PACKAGES_MISSING",
    # Rscript が非ゼロ終了した（stderr 末尾を details に、全文は run の bridge.log）。
    "MSLM_RUN_FAILED",
    # Rscript が制限時間内に終わらなかった（プロセスツリーごと終了済み）。
    "MSLM_TIMEOUT",
    # JSON-LD が無い・壊れている・schemaVersion が想定外（上流の仕様変更の検知）。
    "MSLM_OUTPUT_INVALID",
    # PDF 描画だけが失敗した（ブラウザ不在・unpkg 不達）。pathway_project の結果は有効。
    "PATHWAY_RENDER_FAILED",
    # 引数を直して呼び直すしかない誤り（未知のパスウェイ名・群が小さい・鎖がデータに無い）。
    "PATHWAY_BAD_REQUEST",
})


def pathway_error(code: str, message: str, details: dict | None = None) -> str:
    """パスウェイ投影（MSLipidMapper 連携）固有の確定エラーを返す。"""
    if code not in PATHWAY_ERROR_CODES:
        raise ValueError(f"未知のパスウェイエラーコード: {code!r}")
    if not message:
        raise ValueError("message は必須です。")
    payload: dict = {"code": code, "message": message}
    if details is not None:
        payload["details"] = details
    return json_payload({"error": payload})
```

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_config.py tests/test_user_config.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（バックグラウンド）

```bash
git add metabolomix/core/user_config.py lipidmix.example.toml metabolomix/core/mcp_errors.py tests/test_pathway_config.py
git commit -m "feat(pathway): MSLipidMapper の設定キーとパスウェイ用エラー封筒を足す"
```

---

### Task 2: 中間 mzTab-M（handoff）

**Files:**
- Create: `metabolomix/pathway/__init__.py`（空。docstring 1 行だけ）
- Create: `metabolomix/pathway/handoff.py`
- Test: `tests/test_pathway_handoff.py`

**Interfaces:**
- Produces:
  - `class HandoffError(Exception)`（属性 `code: str`（`PATHWAY_ERROR_CODES` の値）, `message: str`, `details: dict`）
  - `@dataclass HandoffTable: sample_names: list[str]; groups: list[str]; names: list[str]; ontologies: list[str]; matrix: np.ndarray  # (n_features, n_samples)、NaN = 欠損; report: dict`
  - `build_from_arf(*, matrix, sample_names, feature_names, sample_meta, excluded_samples, excluded_spots, catalog: dict[int, dict], props, group_a: str, group_b: str, compared: dict[str, list[str]]) -> HandoffTable`
    （`compared` は `{"group_a": [サンプル名…], "group_b": [...]}`。プール解決は呼び出し側（tools）が既存 `_pool_group_labels` で行う）
  - `build_from_dataset(*, pp_matrix, pp_sample_names, pp_feature_names, roles, sample_names, sample_ids, metadata_rows, annotations: dict[str, dict], group_a: list[str], group_b: list[str], group_a_label: str, group_b_label: str) -> HandoffTable`
  - `write_mztab(table: HandoffTable, path: Path, *, normalization: str) -> None`

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pathway_handoff.py
import math

import numpy as np
import pytest

from metabolomix.pathway import handoff


def _arf_inputs(**over):
    # 4 サンプル（C×2, P×2）＋ QC 1。prop は height と area の 2 種。
    samples = ["C_1", "C_2", "P_1", "P_2", "QC_1"]
    feats = ["Spot_1_height", "Spot_2_height", "Spot_3_height", "Spot_4_height",
             "Spot_5_height", "Spot_1_area"]
    m = np.array([
        # Spot1 Spot2 Spot3 Spot4 Spot5 Spot1area
        [10.0, 1.0, 5.0, 7.0, 3.0, 999.0],
        [12.0, 1.0, 5.0, 7.0, 3.0, 999.0],
        [20.0, 2.0, 6.0, 7.0, np.nan, 999.0],
        [22.0, 2.0, 6.0, 7.0, 4.0, 999.0],
        [15.0, 1.5, 5.5, 7.0, 3.5, 999.0],
    ])
    meta = {s: {"role": "sample", "group": s.split("_")[0]} for s in samples}
    meta["QC_1"] = {"role": "qc", "group": "QC"}
    catalog = {
        1: {"Name": "PC 16:0_18:1", "Ontology": "PC"},
        2: {"Name": "PC 16:0_18:1", "Ontology": "PC"},   # 同名（強度が低い → 落ちる）
        3: {"Name": "Unknown", "Ontology": "Unknown"},
        4: {"Name": "PE 18:0_20:4", "Ontology": "PE"},
        5: {"Name": "Cer 18:1;O2/24:1", "Ontology": "Cer_NS"},
    }
    kw = dict(matrix=m, sample_names=samples, feature_names=feats, sample_meta=meta,
              excluded_samples=set(), excluded_spots=set(), catalog=catalog,
              props=["height", "area"], group_a="C", group_b="P",
              compared={"group_a": ["C_1", "C_2"], "group_b": ["P_1", "P_2"]})
    kw.update(over)
    return kw


def test_arf_keeps_only_biological_samples_with_comparison_labels():
    t = handoff.build_from_arf(**_arf_inputs())
    assert t.sample_names == ["C_1", "C_2", "P_1", "P_2"]
    assert t.groups == ["C", "C", "P", "P"]


def test_arf_pooled_groups_are_labelled_by_the_spec_string():
    kw = _arf_inputs(group_a="ctl", group_b="trt",
                     compared={"group_a": ["C_1", "C_2"], "group_b": ["P_1"]})
    t = handoff.build_from_arf(**kw)
    # 比較外の生体試料 P_2 は自分の群 P のまま残る
    assert dict(zip(t.sample_names, t.groups)) == {
        "C_1": "ctl", "C_2": "ctl", "P_1": "trt", "P_2": "P"}


def test_arf_uses_first_prop_only_and_drops_unknown_and_dedups_by_name():
    t = handoff.build_from_arf(**_arf_inputs())
    assert t.names == ["PC 16:0_18:1", "PE 18:0_20:4", "Cer 18:1;O2/24:1"]
    assert t.ontologies == ["PC", "PE", "Cer_NS"]
    assert t.matrix.shape == (3, 4)
    assert t.matrix[0].tolist() == [10.0, 12.0, 20.0, 22.0]   # Spot1(height) が残る
    assert t.report["dropped_unannotated"] == 1
    assert t.report["dropped_duplicate_names"] == 1
    assert t.report["duplicate_name_examples"] == ["PC 16:0_18:1"]
    assert t.report["prop"] == "height"


def test_arf_excluded_spots_and_samples_are_not_written():
    t = handoff.build_from_arf(**_arf_inputs(excluded_spots={4}, excluded_samples={"P_2"}))
    assert "PE 18:0_20:4" not in t.names
    assert "P_2" not in t.sample_names
    assert t.report["dropped_excluded_spots"] == 1


def test_arf_keeps_nan_as_missing():
    t = handoff.build_from_arf(**_arf_inputs())
    cer = t.names.index("Cer 18:1;O2/24:1")
    assert math.isnan(t.matrix[cer, 2])


def test_arf_too_small_group_is_bad_request():
    kw = _arf_inputs(group_a="ctl", group_b="trt",
                     compared={"group_a": ["C_1"], "group_b": ["P_1", "P_2"]})
    with pytest.raises(handoff.HandoffError) as exc:
        handoff.build_from_arf(**kw)
    assert exc.value.code == "PATHWAY_BAD_REQUEST"


def test_uncompared_sample_whose_group_collides_with_a_label_is_not_written():
    # C_2 は比較外だが自分の群名 "C" が group_a の指定名と同じ。比較群に紛れ込ませると
    # 群 C が 2 例に見えて検定が通ってしまう → 書かない（結果 C=1 で BAD_REQUEST）。
    kw = _arf_inputs(compared={"group_a": ["C_1"], "group_b": ["P_1", "P_2"]})
    with pytest.raises(handoff.HandoffError):
        handoff.build_from_arf(**kw)


def _ds_inputs(**over):
    kw = dict(
        pp_matrix=np.array([[1.0, 50.0], [2.0, 60.0], [3.0, 70.0], [4.0, 80.0], [9.0, 9.0]]),
        pp_sample_names=["a1", "a2", "b1", "b2", "qc"],
        pp_feature_names=["11", "12"],
        roles={"a1": "sample", "a2": "sample", "b1": "sample", "b2": "sample", "qc": "qc"},
        sample_names=["a1", "a2", "b1", "b2", "qc"],
        sample_ids=[],
        metadata_rows=None,
        annotations={"11": {"name": "LPC 18:1", "ontology": "LPC", "ambiguous": False},
                     "12": {"name": None, "ambiguous": True}},
        group_a=["a1", "a2"], group_b=["b1", "b2"],
        group_a_label="ctl", group_b_label="trt")
    kw.update(over)
    return kw


def test_dataset_uses_sml_ontology_and_labels():
    t = handoff.build_from_dataset(**_ds_inputs())
    assert t.sample_names == ["a1", "a2", "b1", "b2"]
    assert t.groups == ["ctl", "ctl", "trt", "trt"]
    assert t.names == ["LPC 18:1"]
    assert t.matrix.tolist() == [[1.0, 2.0, 3.0, 4.0]]
    assert t.report["dropped_unannotated"] == 1


def test_dataset_metadata_groups_add_uncompared_samples():
    names = ["a1", "a2", "b1", "b2", "z1", "x1", "qc"]
    rows = [{"sample_id": n, "group": g, "role": "sample", "include": inc}
            for n, g, inc in [("a1", "A", True), ("a2", "A", True), ("b1", "B", True),
                              ("b2", "B", True), ("z1", "Z", True), ("x1", "X", False)]]
    kw = _ds_inputs(
        pp_matrix=np.arange(14, dtype=float).reshape(7, 2) + 1,
        pp_sample_names=names, sample_names=names,
        roles={**{n: "sample" for n in names}, "qc": "qc"},
        metadata_rows=rows)
    t = handoff.build_from_dataset(**kw)
    # z1 は実験情報シートの群 Z で載る。x1 は include=False、qc は役割で落ちる。
    assert t.sample_names == ["a1", "a2", "b1", "b2", "z1"]
    assert t.groups == ["ctl", "ctl", "trt", "trt", "Z"]


def test_dataset_too_small_group():
    with pytest.raises(handoff.HandoffError):
        handoff.build_from_dataset(**_ds_inputs(group_b=["b1"]))


def test_write_mztab_matches_upstream_minimal_shape(tmp_path):
    t = handoff.build_from_arf(**_arf_inputs())
    out = tmp_path / "input.mztab"
    handoff.write_mztab(t, out, normalization="pqn")
    lines = out.read_text(encoding="utf-8").splitlines()
    mtd = [l.split("\t") for l in lines if l.startswith("MTD")]
    assert ["MTD", "mzTab-version", "2.0.0-M"] in mtd
    assert ["MTD", "comment[1]", "Normalization: pqn"] in mtd
    assert ["MTD", "assay[1]", "C_1"] in mtd
    assert ["MTD", "study_variable[1]", "C"] in mtd
    assert ["MTD", "study_variable[1]-assay_refs", "assay[1]|assay[2]"] in mtd
    assert ["MTD", "study_variable[2]-assay_refs", "assay[3]|assay[4]"] in mtd
    smh = next(l.split("\t") for l in lines if l.startswith("SMH"))
    assert smh[:9] == ["SMH", "SML_ID", "SMF_ID_REFS", "database_identifier",
                       "chemical_formula", "smiles", "inchi", "chemical_name", "uri"]
    assert smh[-1] == "opt_global_Ontology"
    assert smh[14:18] == [f"abundance_assay[{i}]" for i in range(1, 5)]
    sml = [l.split("\t") for l in lines if l.startswith("SML")]
    assert sml[0][7] == "PC 16:0_18:1" and sml[0][-1] == "PC"
    cer = next(r for r in sml if r[7] == "Cer 18:1;O2/24:1")
    assert cer[16] == "null"   # P_1 の NaN
```

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_handoff.py -q`
Expected: FAIL（`ModuleNotFoundError: metabolomix.pathway`）

- [ ] **Step 3: 実装**

```python
# metabolomix/pathway/__init__.py
"""MSLipidMapper を通したパスウェイ投影（spec 2026-10-14）。"""
```

```python
# metabolomix/pathway/handoff.py
"""ARF / DatasetState → 上流 write_mztab_m と同形の最小 mzTab-M（純ロジック）。

群・除外・前処理の決定権はこちらに置く（spec §2）。上流は mzTab-M の
`study_variable` を群に、`opt_global_Ontology` をサブクラスに読む。
session_state / mcp_core を import しない。
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_UNANNOTATED = {"", "unknown", "others", "null"}
_SPOT_RE = re.compile(r"^Spot_(\d+)_(.+)$")
_MIN_GROUP = 2


class HandoffError(Exception):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details or {}


@dataclass
class HandoffTable:
    sample_names: list[str]
    groups: list[str]
    names: list[str]
    ontologies: list[str]
    matrix: np.ndarray            # (n_features, n_samples)。NaN = 欠損
    report: dict = field(default_factory=dict)


def _check_group_sizes(groups: list[str], label_a: str, label_b: str) -> None:
    n_a, n_b = groups.count(label_a), groups.count(label_b)
    if n_a < _MIN_GROUP or n_b < _MIN_GROUP:
        raise HandoffError("PATHWAY_BAD_REQUEST",
                           f"各群 2 例以上が要ります（{label_a}={n_a}, {label_b}={n_b}）。",
                           {"n_a": n_a, "n_b": n_b})


def _annotated(name, ontology) -> bool:
    return (str(name or "").strip().lower() not in _UNANNOTATED
            and str(ontology or "").strip().lower() not in _UNANNOTATED)


def _dedup_by_name(names, ontologies, rows):
    """同名の行は平均強度が最大の 1 行だけ残す（sum の二重計上を避ける。spec §4）。"""
    best: dict[str, int] = {}
    for i, name in enumerate(names):
        mean = float(np.nanmean(rows[i])) if np.isfinite(rows[i]).any() else -math.inf
        j = best.get(name)
        if j is None:
            best[name] = i
            continue
        mean_j = float(np.nanmean(rows[j])) if np.isfinite(rows[j]).any() else -math.inf
        if mean > mean_j:
            best[name] = i
    keep = sorted(best.values())
    dup_names = sorted({n for n in names if names.count(n) > 1})
    return keep, len(names) - len(keep), dup_names[:5]


def _finish(sample_names, groups, names, ontologies, rows, report) -> HandoffTable:
    keep, n_dup, examples = _dedup_by_name(names, ontologies, rows)
    report["dropped_duplicate_names"] = n_dup
    report["duplicate_name_examples"] = examples
    matrix = (np.vstack([rows[i] for i in keep]) if keep
              else np.empty((0, len(sample_names))))
    if not keep:
        raise HandoffError("PATHWAY_BAD_REQUEST",
                           "サブクラス注釈の付いた特徴が 1 件もありません。", dict(report))
    return HandoffTable(sample_names=sample_names, groups=groups,
                        names=[names[i] for i in keep],
                        ontologies=[ontologies[i] for i in keep],
                        matrix=matrix, report=report)


def build_from_arf(*, matrix, sample_names, feature_names, sample_meta,
                   excluded_samples, excluded_spots, catalog, props,
                   group_a: str, group_b: str, compared: dict) -> HandoffTable:
    matrix = np.asarray(matrix, dtype=float)          # (n_samples, n_features)
    prop_list = [props] if isinstance(props, str) else list(props or [])
    prop = prop_list[0] if prop_list else None
    a_set, b_set = set(compared["group_a"]), set(compared["group_b"])

    col_idx, sample_out, groups = [], [], []
    for i, name in enumerate(sample_names):
        entry = sample_meta.get(name) or {}
        if entry.get("role", "sample") != "sample" or name in excluded_samples:
            continue
        if name in a_set:
            label = group_a
        elif name in b_set:
            label = group_b
        else:
            label = entry.get("group")
            # 比較外の試料の群名が比較の指定名と同じなら書かない（比較群に紛れ込む）
            if not label or str(label) in (group_a, group_b):
                continue
        col_idx.append(i)
        sample_out.append(name)
        groups.append(str(label))
    _check_group_sizes(groups, group_a, group_b)

    names, ontologies, rows = [], [], []
    report = {"source": "arf", "prop": prop, "dropped_unannotated": 0,
              "dropped_excluded_spots": 0}
    for j, feat in enumerate(feature_names):
        m = _SPOT_RE.match(str(feat))
        if not m or (prop is not None and m.group(2) != prop):
            continue
        spot_id = int(m.group(1))
        if spot_id in excluded_spots:
            report["dropped_excluded_spots"] += 1
            continue
        spot = catalog.get(spot_id) or {}
        name, onto = spot.get("Name"), spot.get("Ontology")
        if not _annotated(name, onto):
            report["dropped_unannotated"] += 1
            continue
        names.append(str(name).strip())
        ontologies.append(str(onto).strip())
        rows.append(matrix[col_idx, j])
    return _finish(sample_out, groups, names, ontologies, rows, report)


def build_from_dataset(*, pp_matrix, pp_sample_names, pp_feature_names, roles,
                       sample_names, sample_ids, metadata_rows, annotations,
                       group_a, group_b, group_a_label: str,
                       group_b_label: str) -> HandoffTable:
    matrix = np.asarray(pp_matrix, dtype=float)       # (n_samples, n_features)
    a_set, b_set = set(group_a), set(group_b)
    meta_group: dict[str, str] = {}
    if metadata_rows:
        id_of = (dict(zip(sample_names, sample_ids)) if sample_ids
                 else {n: n for n in sample_names})
        by_id = {str(r.get("sample_id")): r for r in metadata_rows}
        for name in sample_names:
            row = by_id.get(str(id_of.get(name, name)))
            if (row and row.get("include", True) and row.get("role", "sample") == "sample"
                    and row.get("group")):
                meta_group[name] = str(row["group"])

    col_idx, sample_out, groups = [], [], []
    for i, name in enumerate(pp_sample_names):
        if roles.get(name, "sample") != "sample":
            continue
        if name in a_set:
            label = group_a_label
        elif name in b_set:
            label = group_b_label
        elif name in meta_group and meta_group[name] not in (group_a_label, group_b_label):
            label = meta_group[name]
        else:
            continue
        col_idx.append(i)
        sample_out.append(name)
        groups.append(label)
    _check_group_sizes(groups, group_a_label, group_b_label)

    names, ontologies, rows = [], [], []
    report = {"source": "dataset", "dropped_unannotated": 0}
    for j, fid in enumerate(pp_feature_names):
        ann = annotations.get(str(fid)) or {}
        if ann.get("ambiguous") or not _annotated(ann.get("name"), ann.get("ontology")):
            report["dropped_unannotated"] += 1
            continue
        names.append(str(ann["name"]).strip())
        ontologies.append(str(ann["ontology"]).strip())
        rows.append(matrix[col_idx, j])
    return _finish(sample_out, groups, names, ontologies, rows, report)


def _cell(v) -> str:
    if v is None or not math.isfinite(float(v)):
        return "null"
    return repr(float(v))


def _clean(s: str) -> str:
    return re.sub(r"[\t\r\n]+", " ", str(s))


def write_mztab(table: HandoffTable, path: Path, *, normalization: str) -> None:
    """上流 R `write_mztab_m` と同じ列構成で書く（R/mztab_m_import.R）。"""
    levels: list[str] = []
    for g in table.groups:
        if g not in levels:
            levels.append(g)
    lines = [["MTD", "mzTab-version", "2.0.0-M"],
             ["MTD", "mzTab-ID", "metabolomix-pathway"],
             ["MTD", "software[1]", "[,,metabolomix,]"],
             ["MTD", "small_molecule-quantification_unit", "[,,abundance,]"],
             ["MTD", "comment[1]", f"Normalization: {_clean(normalization)}"]]
    for i, name in enumerate(table.sample_names, start=1):
        lines.append(["MTD", f"assay[{i}]", _clean(name)])
    for k, level in enumerate(levels, start=1):
        refs = "|".join(f"assay[{i}]" for i, g in enumerate(table.groups, start=1) if g == level)
        lines += [["MTD", f"study_variable[{k}]", _clean(level)],
                  ["MTD", f"study_variable[{k}]-assay_refs", refs],
                  ["MTD", f"study_variable[{k}]-description", _clean(level)]]
    n = len(table.sample_names)
    lines.append(["SMH", "SML_ID", "SMF_ID_REFS", "database_identifier",
                  "chemical_formula", "smiles", "inchi", "chemical_name", "uri",
                  "theoretical_neutral_mass", "adduct_ions", "reliability",
                  "best_id_confidence_measure", "best_id_confidence_value",
                  *[f"abundance_assay[{i}]" for i in range(1, n + 1)],
                  "opt_global_Ontology"])
    for r, (name, onto) in enumerate(zip(table.names, table.ontologies), start=1):
        lines.append(["SML", str(r), str(r), "null", "null", "null", "null",
                      _clean(name), "null", "null", "null", "null", "null", "null",
                      *[_cell(v) for v in table.matrix[r - 1]], _clean(onto)])
    path.write_text("\n".join("\t".join(l) for l in lines) + "\n", encoding="utf-8")
```

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_handoff.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（バックグラウンド）

```bash
git add metabolomix/pathway/__init__.py metabolomix/pathway/handoff.py tests/test_pathway_handoff.py
git commit -m "feat(pathway): ARF / DatasetState から上流と同形の最小 mzTab-M を書く"
```

---

### Task 3: DatasetState が SML の `opt_global_Ontology` を持つ

**Files:**
- Modify: `metabolomix/mztab/dataset_state.py`（`_SML_TEXT_FIELDS`）
- Modify: `docs/output_format/mztab.md`（`feature_annotations` の項に `ontology` を一行追加）
- Test: `tests/test_dataset_state.py`（末尾に追加）

**Interfaces:**
- Produces: `ds.feature_annotations[fid]["ontology"]: str | None`（Task 8 の dataset 経路が読む）。

- [ ] **Step 1: 失敗するテストを書く**（`tests/test_dataset_state.py` 末尾。既存の `_build_text` と `_SML_ONLY_MZTAB` を流用し、SMH/SML に列を足した派生文字列を作る）

```python
def test_sml_opt_global_ontology_is_kept_as_ontology(tmp_path):
    text = _SML_ONLY_MZTAB.replace(
        "best_id_confidence_value", "best_id_confidence_value\topt_global_Ontology", 1)
    # SML 行の末尾にサブクラスを足す（SMH に列を足したので各 SML 行にも 1 列）。
    rows = []
    for line in text.splitlines():
        if line.startswith("SML\t"):
            line = line + ("\tPC" if line.split("\t")[1] == "1" else "\tnull")
        rows.append(line)
    ds = _build_text(tmp_path, "\n".join(rows) + "\n", name="onto.mzTab")
    assert ds.feature_annotations["1"]["ontology"] == "PC"
```

実装前に `_SML_ONLY_MZTAB` の SMH 行を開き、`best_id_confidence_value` が SMH の最終列か確かめる。最終列でなければ `replace` の対象を SMH の最終列名に変える（列の並びは SML 行の末尾追記と対応していればよい）。`_build_text` の引数名も実物に合わせる。

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_dataset_state.py -q -k ontology`
Expected: FAIL（`KeyError: 'ontology'`）

- [ ] **Step 3: 実装**

`_SML_TEXT_FIELDS` に 1 行追加:

```python
    # MS-DIAL の SML 拡張列。脂質サブクラス（パスウェイ投影の照合キー）。
    "ontology": "opt_global_Ontology",
```

`docs/output_format/mztab.md` の `feature_annotations` を説明している箇所に追記:

```markdown
- `ontology` — SML の `opt_global_Ontology`（MS-DIAL の脂質サブクラス）。列が無いファイルでは `None`。
  `pathway_project` の照合キーになる。差次的エクスポートの列には出ない。
```

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_dataset_state.py tests/test_mztab_tools.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（バックグラウンド）

```bash
git add metabolomix/mztab/dataset_state.py docs/output_format/mztab.md tests/test_dataset_state.py
git commit -m "feat(mztab): SML の opt_global_Ontology を feature_annotations の ontology として持つ"
```

---

### Task 4: JSON-LD の読み取り

**Files:**
- Create: `metabolomix/pathway/jsonld.py`
- Test: `tests/test_pathway_jsonld.py`

**Interfaces:**
- Produces:
  - `class JsonldError(Exception)`（`code="MSLM_OUTPUT_INVALID"`, `message`）
  - `SCHEMA_VERSION = "2.0.0"`
  - `@dataclass Scope: id: str; node_id: str; scope_type: str  # "subclass" | "acyl_subset"; chains: tuple[str, ...]; n_features: int`
  - `@dataclass PathwayGraph: node_labels: dict[str, str]; node_kinds: dict[str, str]  # "lipid" | "entity" | "junction"; edges: list[tuple[str, str]]; scopes: dict[str, Scope]; sample_groups: dict[str, str]  # identifier -> group; values: dict[str, dict[str, float | None]]  # scope_id -> identifier -> value; unmatched_subclasses: list[str]`
  - `load_graph(path: Path) -> PathwayGraph`
  - `adjacent_lipid_pairs(graph: PathwayGraph) -> list[tuple[str, str, str]]`（`(node_id_1, node_id_2, via)`、`via` は `"direct"` か `"junction"`。同じ組は 1 回、`node_id_1 < node_id_2`）

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pathway_jsonld.py
import json

import pytest

from metabolomix.pathway import jsonld

B = "https://x/d/"
V = "https://systemomicslab.github.io/MSLipidMapper/ontology/"


def _ref(i):
    return {"@id": B + i}


def _graph(extra=(), schema="2.0.0"):
    g = [
        {"@id": B + "ds", "@type": "mslm:Dataset", "schemaVersion": schema,
         "identifier": "t", "pathway": _ref("pw")},
        {"@id": B + "pw", "@type": "mslm:Pathway"},
        {"@id": B + "sub/PC", "@type": "mslm:MSDialMappableLipidEntity", "label": "PC"},
        {"@id": B + "sub/PE", "@type": "mslm:MSDialMappableLipidEntity", "label": "PE"},
        {"@id": B + "sub/LPC", "@type": "mslm:MSDialMappableLipidEntity", "label": "LPC"},
        {"@id": B + "sub/TG", "@type": "mslm:MSDialMappableLipidEntity", "label": "TG"},
        {"@id": B + "n/PC", "@type": "mslm:PathwayLipidNode", "label": "PC", "representsEntity": _ref("sub/PC")},
        {"@id": B + "n/PE", "@type": "mslm:PathwayLipidNode", "label": "PE", "representsEntity": _ref("sub/PE")},
        {"@id": B + "n/LPC", "@type": "mslm:PathwayLipidNode", "label": "LPC", "representsEntity": _ref("sub/LPC")},
        {"@id": B + "n/PA", "@type": "mslm:PathwayEntity", "label": "PA"},
        {"@id": B + "n/DG", "@type": "mslm:PathwayLipidNode", "label": "DG"},
        {"@id": B + "n/j1", "@type": "mslm:PathwayJunction"},
        # PC -> j1 -> PE, j1 -> LPC（分岐）／ PE -> PA -> DG（PA で隣接が切れる）／ PC -> LPC 直結
        {"@id": B + "c1", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PC"), "targetEntity": _ref("n/j1")},
        {"@id": B + "c2", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/j1"), "targetEntity": _ref("n/PE")},
        {"@id": B + "c3", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/j1"), "targetEntity": _ref("n/LPC")},
        {"@id": B + "c4", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PE"), "targetEntity": _ref("n/PA")},
        {"@id": B + "c5", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PA"), "targetEntity": _ref("n/DG")},
        {"@id": B + "c6", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PC"), "targetEntity": _ref("n/LPC")},
        {"@id": B + "s/1", "@type": "mslm:Sample", "identifier": "S1", "sampleGroup": "C"},
        {"@id": B + "s/2", "@type": "mslm:Sample", "identifier": "S2", "sampleGroup": "P"},
        {"@id": B + "n/PC/scope/whole", "@type": "mslm:LipidMeasurementScope", "onPathwayNode": _ref("n/PC"),
         "scopeType": "subclass", "acylChainSignature": [], "memberFeature": [_ref("f/1"), _ref("f/2")]},
        {"@id": B + "n/PC/scope/acyl-18%3A0", "@type": "mslm:LipidMeasurementScope", "onPathwayNode": _ref("n/PC"),
         "scopeType": "acyl_subset", "acylChainSignature": ["18:0"], "memberFeature": [_ref("f/2")]},
        {"@id": B + "m1", "@type": "mslm:SampleMeasurement", "measuredEntity": _ref("n/PC/scope/whole"),
         "inSample": _ref("s/1"), "value": 15, "valueStatus": {"@id": V + "AvailableValue"}},
        {"@id": B + "m2", "@type": "mslm:SampleMeasurement", "measuredEntity": _ref("n/PC/scope/whole"),
         "inSample": _ref("s/2"), "value": None, "valueStatus": {"@id": V + "MissingValue"}},
        {"@id": B + "m3", "@type": "mslm:SampleMeasurement", "measuredEntity": _ref("f/1"),
         "inSample": _ref("s/1"), "value": 10},
        *extra,
    ]
    return {"@context": {}, "@graph": g}


def _write(tmp_path, doc):
    p = tmp_path / "remodeling-measurements.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def test_nodes_scopes_samples_and_values(tmp_path):
    g = jsonld.load_graph(_write(tmp_path, _graph()))
    assert g.node_kinds[B + "n/PC"] == "lipid"
    assert g.node_kinds[B + "n/PA"] == "entity"
    assert g.node_kinds[B + "n/j1"] == "junction"
    whole = g.scopes[B + "n/PC/scope/whole"]
    assert (whole.scope_type, whole.chains, whole.n_features) == ("subclass", (), 2)
    assert g.scopes[B + "n/PC/scope/acyl-18%3A0"].chains == ("18:0",)
    assert g.sample_groups == {"S1": "C", "S2": "P"}
    assert g.values[B + "n/PC/scope/whole"] == {"S1": 15.0, "S2": None}
    # 分子種（feature）の測定値はスコープ値に混ぜない
    assert B + "f/1" not in g.values


def test_unmatched_subclasses(tmp_path):
    g = jsonld.load_graph(_write(tmp_path, _graph()))
    assert g.unmatched_subclasses == ["TG"]


def test_adjacency_through_junctions_and_not_through_entities(tmp_path):
    g = jsonld.load_graph(_write(tmp_path, _graph()))
    pairs = {(g.node_labels[a], g.node_labels[b], via)
             for a, b, via in jsonld.adjacent_lipid_pairs(g)}
    assert ("LPC", "PC", "direct") in pairs or ("PC", "LPC", "direct") in pairs
    labels = {frozenset((a, b)) for a, b, _ in pairs}
    assert frozenset(("PC", "PE")) in labels
    assert frozenset(("LPC", "PE")) in labels      # 分岐した接合ノードの枝同士
    assert frozenset(("PE", "DG")) not in labels   # PA（entity）で切れる
    # PC–LPC は direct と junction の両方で結ばれているが 1 回だけ
    assert sum(1 for a, b, _ in pairs if {a, b} == {"PC", "LPC"}) == 1


def test_wrong_schema_version_is_rejected(tmp_path):
    with pytest.raises(jsonld.JsonldError) as exc:
        jsonld.load_graph(_write(tmp_path, _graph(schema="3.0.0")))
    assert exc.value.code == "MSLM_OUTPUT_INVALID"


def test_broken_json_is_rejected(tmp_path):
    p = tmp_path / "x.json"
    p.write_text("{", encoding="utf-8")
    with pytest.raises(jsonld.JsonldError):
        jsonld.load_graph(p)
```

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_jsonld.py -q`
Expected: FAIL（`ImportError`）

- [ ] **Step 3: 実装**

```python
# metabolomix/pathway/jsonld.py
"""MSLipidMapper の pathway measurements JSON-LD（schema 2.0.0）を読む（純ロジック）。

上流文書: docs/pathway-measurements.md（MSLipidMapper リポジトリ）。型は `mslm:` 接頭辞付きの
compact 形で来る。参照は {"@id": ...}。値の欠損は JSON null（MissingValue）。
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA_VERSION = "2.0.0"
_KIND = {"mslm:PathwayLipidNode": "lipid", "mslm:PathwayEntity": "entity",
         "mslm:PathwayJunction": "junction"}


class JsonldError(Exception):
    code = "MSLM_OUTPUT_INVALID"

    def __init__(self, message: str):
        super().__init__(message)
        self.message = message


@dataclass
class Scope:
    id: str
    node_id: str
    scope_type: str
    chains: tuple[str, ...]
    n_features: int


@dataclass
class PathwayGraph:
    node_labels: dict[str, str] = field(default_factory=dict)
    node_kinds: dict[str, str] = field(default_factory=dict)
    edges: list[tuple[str, str]] = field(default_factory=list)
    scopes: dict[str, Scope] = field(default_factory=dict)
    sample_groups: dict[str, str] = field(default_factory=dict)
    values: dict[str, dict[str, float | None]] = field(default_factory=dict)
    unmatched_subclasses: list[str] = field(default_factory=list)


def _ref(obj, key) -> str | None:
    value = obj.get(key)
    return value.get("@id") if isinstance(value, dict) else None


def load_graph(path: Path) -> PathwayGraph:
    try:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        items = doc["@graph"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise JsonldError(f"JSON-LD を読めません: {path} ({exc})") from exc
    datasets = [x for x in items if x.get("@type") == "mslm:Dataset"]
    version = datasets[0].get("schemaVersion") if datasets else None
    if version != SCHEMA_VERSION:
        raise JsonldError(f"schemaVersion が {version!r} です（想定 {SCHEMA_VERSION}）。"
                          "上流 MSLipidMapper の出力仕様が変わった可能性があります。")

    g = PathwayGraph()
    sample_uri_to_id: dict[str, str] = {}
    subclass_labels: dict[str, str] = {}
    represented: set[str] = set()
    measurements = []
    for x in items:
        t = x.get("@type")
        if t in _KIND:
            g.node_kinds[x["@id"]] = _KIND[t]
            g.node_labels[x["@id"]] = str(x.get("label") or "")
            rep = _ref(x, "representsEntity")
            if rep:
                represented.add(rep)
        elif t == "mslm:PathwayConnection":
            src, tgt = _ref(x, "sourceEntity"), _ref(x, "targetEntity")
            if src and tgt:
                g.edges.append((src, tgt))
        elif t == "mslm:LipidMeasurementScope":
            g.scopes[x["@id"]] = Scope(
                id=x["@id"], node_id=_ref(x, "onPathwayNode") or "",
                scope_type=str(x.get("scopeType") or ""),
                chains=tuple(x.get("acylChainSignature") or ()),
                n_features=len(x.get("memberFeature") or ()))
        elif t == "mslm:Sample":
            sample_uri_to_id[x["@id"]] = str(x.get("identifier"))
            g.sample_groups[str(x.get("identifier"))] = str(x.get("sampleGroup") or "")
        elif t == "mslm:MSDialMappableLipidEntity":
            subclass_labels[x["@id"]] = str(x.get("label") or "")
        elif t == "mslm:SampleMeasurement":
            measurements.append(x)

    values: dict[str, dict[str, float | None]] = defaultdict(dict)
    for x in measurements:
        entity = _ref(x, "measuredEntity")
        if entity not in g.scopes:
            continue                       # 分子種（feature）の値はスコープ値に混ぜない
        sample = sample_uri_to_id.get(_ref(x, "inSample") or "")
        if sample is None:
            continue
        v = x.get("value")
        values[entity][sample] = None if v is None else float(v)
    g.values = dict(values)
    g.unmatched_subclasses = sorted(label for uri, label in subclass_labels.items()
                                    if uri not in represented)
    return g


def adjacent_lipid_pairs(graph: PathwayGraph) -> list[tuple[str, str, str]]:
    """接合ノードだけを挟んでつながる脂質ノードの組（向きは問わない）。"""
    neighbors: dict[str, set[str]] = defaultdict(set)
    direct: set[frozenset] = set()
    for a, b in graph.edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
        if graph.node_kinds.get(a) == "lipid" and graph.node_kinds.get(b) == "lipid":
            direct.add(frozenset((a, b)))
    via_junction: set[frozenset] = set()
    for start, kind in graph.node_kinds.items():
        if kind != "lipid":
            continue
        seen, stack = {start}, [n for n in neighbors[start]
                                if graph.node_kinds.get(n) == "junction"]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            for nxt in neighbors[node]:
                k = graph.node_kinds.get(nxt)
                if k == "junction" and nxt not in seen:
                    stack.append(nxt)
                elif k == "lipid" and nxt != start:
                    via_junction.add(frozenset((start, nxt)))
    out = []
    for pair in direct | via_junction:
        a, b = sorted(pair)
        out.append((a, b, "direct" if pair in direct else "junction"))
    return sorted(out)
```

（テストの `("LPC","PC","direct") in pairs or ...` はラベル順に依存しないための書き方。`direct` が優先されることを見ている。）

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_jsonld.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（バックグラウンド）

```bash
git add metabolomix/pathway/jsonld.py tests/test_pathway_jsonld.py
git commit -m "feat(pathway): MSLipidMapper の JSON-LD を読み、接合ノード経由の隣接を求める"
```

---

### Task 5: 2 群比較と逆向き隣接ペア

**Files:**
- Create: `metabolomix/pathway/compare.py`
- Test: `tests/test_pathway_compare.py`

**Interfaces:**
- Consumes: `jsonld.PathwayGraph`, `jsonld.adjacent_lipid_pairs`, `metabolomix.analysis.differential.two_group_test` / `bh_fdr`。
- Produces:
  - `compare_scopes(graph, group_a: str, group_b: str, q_threshold: float) -> list[dict]`
    各 dict のキー: `node scope node_id n_features mean_a mean_b log2fc q call`（`scope` は `"whole"` か鎖 `"18:0"`、`call` は `"up" | "down" | "ns" | "NA"`）。並びはスコープ種別（whole が先）→ ノード名。
  - `opposite_pairs(graph, rows: list[dict]) -> tuple[list[dict], int]`
    （`whole` スコープだけで判定。戻り値 `(pairs, n_undeterminable)`、pair のキー `node_1 call_1 log2fc_1 node_2 call_2 log2fc_2 via`）
  - `rows_to_tsv(rows: list[dict], columns: list[str]) -> str`（float は有効 4 桁）

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pathway_compare.py
from metabolomix.pathway import compare
from metabolomix.pathway.jsonld import PathwayGraph, Scope


def _graph(values_by_node: dict[str, list], edges, kinds=None):
    """values_by_node: label -> [S1..S6]（S1-3=C, S4-6=P）。None は欠損。"""
    g = PathwayGraph()
    for label in values_by_node:
        nid = f"n/{label}"
        g.node_kinds[nid] = "lipid"
        g.node_labels[nid] = label
        sid = f"{nid}/scope/whole"
        g.scopes[sid] = Scope(sid, nid, "subclass", (), 3)
        g.values[sid] = {f"S{i}": v for i, v in enumerate(values_by_node[label], start=1)}
    for nid, kind in (kinds or {}).items():
        g.node_kinds[nid] = kind
        g.node_labels[nid] = ""
    g.edges = [(f"n/{a}" if not a.startswith("j") else a,
                f"n/{b}" if not b.startswith("j") else b) for a, b in edges]
    g.sample_groups = {f"S{i}": ("C" if i <= 3 else "P") for i in range(1, 7)}
    return g


UP = [10, 11, 10.5, 100, 110, 105]
DOWN = [100, 110, 105, 10, 11, 10.5]
FLAT = [10, 50, 30, 12, 48, 31]


def test_up_down_ns_calls_and_direction():
    g = _graph({"PC": UP, "LPC": DOWN, "PE": FLAT}, [("PC", "LPC")])
    rows = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert rows["PC"]["call"] == "up" and rows["PC"]["log2fc"] > 0
    assert rows["LPC"]["call"] == "down"
    assert rows["PE"]["call"] == "ns"
    assert rows["PC"]["scope"] == "whole" and rows["PC"]["n_features"] == 3


def test_all_missing_in_one_group_is_na_and_not_in_bh():
    g = _graph({"PC": UP, "TG": [1, 2, 3, None, None, None]}, [])
    rows = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert rows["TG"]["call"] == "NA"
    assert rows["TG"]["q"] is None
    # TG が BH の分母に入らない＝ PC の q は単独検定の p と一致
    single = compare.compare_scopes(_graph({"PC": UP}, []), "C", "P", 0.05)[0]
    assert rows["PC"]["q"] == single["q"]


def test_bh_is_per_scope_kind():
    g = _graph({"PC": UP}, [])
    g.scopes["n/PC/scope/acyl-18:0"] = Scope("n/PC/scope/acyl-18:0", "n/PC", "acyl_subset", ("18:0",), 1)
    g.values["n/PC/scope/acyl-18:0"] = {f"S{i}": v for i, v in enumerate(FLAT, start=1)}
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    assert [r["scope"] for r in rows] == ["whole", "18:0"]
    whole_alone = compare.compare_scopes(_graph({"PC": UP}, []), "C", "P", 0.05)[0]
    assert rows[0]["q"] == whole_alone["q"]   # 鎖スコープが whole の BH に混ざらない


def test_opposite_pairs_needs_both_significant():
    g = _graph({"PC": UP, "LPC": DOWN, "PE": FLAT, "PS": UP},
               [("PC", "j1"), ("j1", "LPC"), ("PC", "PE"), ("PS", "PC")],
               kinds={"j1": "junction"})
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    pairs, undeterminable = compare.opposite_pairs(g, rows)
    assert [(p["node_1"], p["node_2"], p["via"]) for p in pairs] == [("LPC", "PC", "junction")]
    assert undeterminable == 1          # PC–PE（PE が ns）。PS–PC は同じ向きで不掲載・判定不能でもない


def test_tsv_rounds_floats():
    tsv = compare.rows_to_tsv([{"a": 1.23456789, "b": None, "c": "x"}], ["a", "b", "c"])
    assert tsv == "a\tb\tc\n1.235\t\tx"
```

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_compare.py -q`
Expected: FAIL（`ImportError`）

- [ ] **Step 3: 実装**

```python
# metabolomix/pathway/compare.py
"""JSON-LD のスコープ値で 2 群比較し、隣接ノードの逆向きを拾う（純ロジック）。

検定は既存の v1 経路と同じ: log2(x+1) 空間の Welch（`differential.two_group_test`）。
BH はパスウェイ（＝1 グラフ）× スコープ種別（whole / 鎖ごと）の単位でかける。
"""
from __future__ import annotations

import math

import numpy as np

from metabolomix.analysis import differential
from metabolomix.pathway.jsonld import PathwayGraph, adjacent_lipid_pairs


def _scope_key(scope) -> str:
    return "whole" if scope.scope_type == "subclass" else "+".join(scope.chains)


def _num(v):
    return None if v is None or not math.isfinite(float(v)) else float(v)


def compare_scopes(graph: PathwayGraph, group_a: str, group_b: str,
                   q_threshold: float) -> list[dict]:
    samples = [s for s, g in graph.sample_groups.items() if g in (group_a, group_b)]
    labels = [graph.sample_groups[s] for s in samples]
    by_kind: dict[str, list] = {}
    for scope in graph.scopes.values():
        by_kind.setdefault(_scope_key(scope), []).append(scope)

    out: list[dict] = []
    kinds = ["whole"] + sorted(k for k in by_kind if k != "whole")
    for kind in kinds:
        scopes = sorted(by_kind.get(kind, []), key=lambda s: graph.node_labels.get(s.node_id, ""))
        if not scopes:
            continue
        matrix = np.array([[np.nan if graph.values.get(s.id, {}).get(name) is None
                            else graph.values[s.id][name] for s in scopes]
                           for name in samples], dtype=float)
        results = differential.two_group_test(matrix, [s.id for s in scopes], labels,
                                              group_a, group_b, log_transform=True)
        qs = differential.bh_fdr([r["p"] for r in results])
        for scope, r, q in zip(scopes, results, qs):
            q = _num(q)
            fc = _num(r["log2fc"])
            if q is None or fc is None:
                call = "NA"
            elif q < q_threshold:
                call = "up" if fc > 0 else "down"
            else:
                call = "ns"
            out.append({"node": graph.node_labels.get(scope.node_id, ""),
                        "scope": kind, "node_id": scope.node_id,
                        "n_features": scope.n_features,
                        "mean_a": _num(r["mean_a"]), "mean_b": _num(r["mean_b"]),
                        "log2fc": fc, "q": q, "call": call})
    return out


def opposite_pairs(graph: PathwayGraph, rows: list[dict]) -> tuple[list[dict], int]:
    whole = {r["node_id"]: r for r in rows if r["scope"] == "whole"}
    pairs, undeterminable = [], 0
    for a, b, via in adjacent_lipid_pairs(graph):
        ra, rb = whole.get(a), whole.get(b)
        if ra is None or rb is None:
            continue                     # 値の付かないノード同士は比較の対象外
        calls = {ra["call"], rb["call"]}
        if calls & {"ns", "NA"}:
            undeterminable += 1
            continue
        if ra["call"] != rb["call"]:
            first, second = sorted((ra, rb), key=lambda r: r["node"])
            pairs.append({"node_1": first["node"], "call_1": first["call"],
                          "log2fc_1": first["log2fc"],
                          "node_2": second["node"], "call_2": second["call"],
                          "log2fc_2": second["log2fc"], "via": via})
    pairs.sort(key=lambda p: (p["node_1"], p["node_2"]))
    return pairs, undeterminable


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def rows_to_tsv(rows: list[dict], columns: list[str]) -> str:
    lines = ["\t".join(columns)]
    lines += ["\t".join(_fmt(r.get(c)) for c in columns) for r in rows]
    return "\n".join(lines)
```

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_compare.py -q`
Expected: PASS。`test_opposite_pairs_needs_both_significant` の PS–PC は `PS` と `PC` がどちらも `up` で、掲載されず判定不能にも数えないこと。

- [ ] **Step 5: コミット**（バックグラウンド）

```bash
git add metabolomix/pathway/compare.py tests/test_pathway_compare.py
git commit -m "feat(pathway): スコープ値の 2 群比較と、逆向きに動く隣接ペアの抽出"
```

---

### Task 6: R ブリッジ（Python ランナー＋ `mslm_bridge.R`）

**Files:**
- Create: `metabolomix/pathway/mslm_bridge.R`
- Create: `metabolomix/pathway/bridge.py`
- Test: `tests/test_pathway_bridge.py`
- Test fixture（テスト内で生成）: 偽 Rscript は Python スクリプト

**Interfaces:**
- Produces:
  - `BRIDGE_SCRIPT: Path`（= `Path(__file__).with_name("mslm_bridge.R")`。テストが差し替える）
  - `PATHWAY_NETWORKS: dict[str, str] = {"remodeling": "remodeling_lipidonly.cyjs", "ceramide": "ceramidepathway_lipidonly.cyjs", "global": "global.cyjs"}`
  - `class BridgeError(Exception)`（`code`, `message`, `details`）
  - `resolve_rscript() -> str`（未設定・不在は `BridgeError("MSLM_NOT_CONFIGURED")`）
  - `upstream_commit() -> str | None`（`[mslipidmapper] repo` があれば `git -C <repo> rev-parse HEAD`、無ければ None）
  - `run_bridge(request: dict, run_dir: Path, *, timeout_s: float = 300.0) -> dict`
    （`run_dir/request.json` を書き、`[rscript, BRIDGE_SCRIPT, request.json]` を `start_owned_process` で起動。結果 `run_dir/bridge-result.json` を dict で返す）

**ブリッジの request / result 契約**（R と Python の両方がこれを守る）:

```
request.json
{"mode": "check" | "measure" | "render",
 "mztab": "<abs path>", "out_dir": "<abs path>",
 "pathways": ["remodeling", ...],
 "networks": {"remodeling": "remodeling_lipidonly.cyjs", ...},
 "acyl_chains": ["18:0"], "aggregation": "sum",
 "normalization": "<recipe text>", "unit": "arbitrary", "dataset_id": "<run_id>",
 "browser": "<path or null>", "plot_type": "dot", "x_order": ["C", "P"],
 "timeout_seconds": 90}

bridge-result.json
{"status": "ok" | "error",
 "code": null | "MSLM_R_PACKAGES_MISSING" | "PATHWAY_RENDER_FAILED" | "PATHWAY_BAD_REQUEST" | "MSLM_RUN_FAILED",
 "message": "...",
 "package_version": "0.1.0.9000",
 "outputs": {"remodeling": "<abs path of json or pdf>", ...},
 "mapped_nodes": {"remodeling": 12, ...}}
```

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pathway_bridge.py
import json
import sys
import textwrap

import pytest

from metabolomix.pathway import bridge

pytestmark = pytest.mark.skipif(sys.platform != "win32", reason="start_owned_process は Windows 専用")


def _fake(tmp_path, body: str):
    """偽ブリッジ: Python スクリプト。argv[1] が request.json。"""
    script = tmp_path / "fake_bridge.py"
    script.write_text(textwrap.dedent("""
        import json, sys, time, pathlib
        req = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
        out = pathlib.Path(req["out_dir"])
    """) + textwrap.dedent(body), encoding="utf-8")
    return script


@pytest.fixture
def configured(monkeypatch, tmp_path):
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", sys.executable)
    return tmp_path


def test_not_configured(monkeypatch):
    monkeypatch.delenv("MSLIPIDMAPPER_RSCRIPT", raising=False)
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.resolve_rscript()
    assert exc.value.code == "MSLM_NOT_CONFIGURED"


def test_configured_but_missing_file(monkeypatch, tmp_path):
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", str(tmp_path / "nope.exe"))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.resolve_rscript()
    assert exc.value.code == "MSLM_NOT_CONFIGURED"


def test_ok_result_is_returned(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        (out / "bridge-result.json").write_text(json.dumps(
            {"status": "ok", "code": None, "message": "", "outputs": {"remodeling": "x.json"},
             "mapped_nodes": {"remodeling": 3}, "package_version": "0.1"}), encoding="utf-8")
    """))
    res = bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert res["mapped_nodes"] == {"remodeling": 3}
    assert json.loads((run_dir / "request.json").read_text(encoding="utf-8"))["mode"] == "measure"


def test_error_result_carries_code(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        (out / "bridge-result.json").write_text(json.dumps(
            {"status": "error", "code": "MSLM_R_PACKAGES_MISSING",
             "message": "MSLipidMapper が入っていません"}), encoding="utf-8")
        sys.exit(2)
    """))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "check", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_R_PACKAGES_MISSING"


def test_nonzero_without_result_is_run_failed_with_log_tail(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        for i in range(40):
            print(f"line {i}", file=sys.stderr)
        sys.exit(1)
    """))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_RUN_FAILED"
    tail = exc.value.details["log_tail"]
    assert tail.splitlines()[-1] == "line 39" and len(tail.splitlines()) == 20
    assert (run_dir / "bridge.log").is_file()


def test_timeout_kills_tree(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, """
        import subprocess
        subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
        time.sleep(60)
    """))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=2)
    assert exc.value.code == "MSLM_TIMEOUT"


def test_zero_exit_without_result_is_output_invalid(configured, monkeypatch):
    run_dir = configured / "run"
    run_dir.mkdir()
    monkeypatch.setattr(bridge, "BRIDGE_SCRIPT", _fake(configured, "pass\n"))
    with pytest.raises(bridge.BridgeError) as exc:
        bridge.run_bridge({"mode": "measure", "out_dir": str(run_dir)}, run_dir, timeout_s=60)
    assert exc.value.code == "MSLM_OUTPUT_INVALID"
```

`test_timeout_kills_tree` は孫プロセスが残らないことを `OwnedProcess.terminate_tree` に任せている（Job Object の KILL_ON_JOB_CLOSE）。孫の生存確認まではしない（既存 `tests/test_process_control.py` が持つ）。

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_bridge.py -q`
Expected: FAIL（`ImportError`）

- [ ] **Step 3: Python 側を実装**

```python
# metabolomix/pathway/bridge.py
"""MSLipidMapper（R パッケージ）を Rscript サブプロセスで呼ぶ。

契約は request.json → bridge-result.json（R 側は mslm_bridge.R）。Console と同じく
自前の Job Object で起動し、timeout ではプロセスツリーごと終了する。再試行しない。
session_state / mcp_core を import しない。
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from metabolomix.core import user_config
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.process_control import start_owned_process

BRIDGE_SCRIPT = Path(__file__).with_name("mslm_bridge.R")
PATHWAY_NETWORKS = {"remodeling": "remodeling_lipidonly.cyjs",
                    "ceramide": "ceramidepathway_lipidonly.cyjs",
                    "global": "global.cyjs"}
_LOG_TAIL_LINES = 20


class BridgeError(Exception):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details or {}


def resolve_rscript() -> str:
    key = "mslipidmapper.rscript"
    setting = user_config.get_setting(key)
    if setting is None:
        raise BridgeError("MSLM_NOT_CONFIGURED",
                          user_config.missing_hint(key, "Rscript（MSLipidMapper 入りの R）"),
                          user_config.describe_missing(key))
    if not Path(setting.value).is_file():
        raise BridgeError("MSLM_NOT_CONFIGURED",
                          f"{user_config.setting_label(setting)} が指す Rscript がありません。",
                          user_config.describe_missing(key, setting))
    return setting.value


def browser_path() -> str | None:
    setting = user_config.get_setting("mslipidmapper.browser")
    return setting.value if setting else None


def upstream_commit() -> str | None:
    setting = user_config.get_setting("mslipidmapper.repo")
    if setting is None:
        return None
    try:
        out = subprocess.run(["git", "-C", setting.value, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def _log_tail(log_path: Path) -> str:
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-_LOG_TAIL_LINES:])


def run_bridge(request: dict, run_dir: Path, *, timeout_s: float = 300.0) -> dict:
    rscript = resolve_rscript()
    run_dir = Path(run_dir)
    request_path = run_dir / "request.json"
    result_path = run_dir / "bridge-result.json"
    log_path = run_dir / "bridge.log"
    result_path.unlink(missing_ok=True)
    request_path.write_text(json.dumps(request, ensure_ascii=False), encoding="utf-8")

    proc = start_owned_process([rscript, str(BRIDGE_SCRIPT), str(request_path)],
                               cwd=run_dir, log_path=log_path)
    try:
        try:
            exit_code = proc.wait(timeout=timeout_s)
        except DomainError as exc:
            if exc.code != "PROCESS_WAIT_TIMEOUT":
                raise
            proc.terminate_tree()
            raise BridgeError("MSLM_TIMEOUT",
                              f"MSLipidMapper が {timeout_s:.0f} 秒以内に終わりませんでした。",
                              {"timeout_s": timeout_s, "log": str(log_path)}) from exc
    finally:
        proc.close()

    result = None
    if result_path.is_file():
        try:
            result = json.loads(result_path.read_text(encoding="utf-8"))
        except ValueError:
            result = None
    if result is not None and result.get("status") == "error":
        raise BridgeError(result.get("code") or "MSLM_RUN_FAILED",
                          result.get("message") or "MSLipidMapper が失敗しました。",
                          {"log_tail": _log_tail(log_path), "log": str(log_path)})
    if exit_code != 0:
        raise BridgeError("MSLM_RUN_FAILED",
                          f"Rscript が終了コード {exit_code} で終わりました。",
                          {"exit_code": exit_code, "log_tail": _log_tail(log_path),
                           "log": str(log_path)})
    if result is None:
        raise BridgeError("MSLM_OUTPUT_INVALID",
                          "Rscript は正常終了しましたが bridge-result.json がありません。",
                          {"log_tail": _log_tail(log_path), "log": str(log_path)})
    return result
```

`DomainError` は `metabolomix/core/atomic_io.py` の定義（`process_control.py` の `OwnedProcess.wait` が timeout で送出する `PROCESS_WAIT_TIMEOUT` もこれ。`exc.code` で判定する——実装前に `DomainError` の属性名が `code` であることを確認）。`start_owned_process` は stdout/stderr を `log_path` へ落とす（`_open_child_stdio`）。

- [ ] **Step 4: R 側を実装**（R はテストで走らない。Task 11 の実データ検証で走る）

```r
# metabolomix/pathway/mslm_bridge.R
# MSLipidMapper（R パッケージ）を呼ぶ薄いブリッジ。契約は bridge.py の docstring と計画書を参照。
# 引数: request.json のパス。結果は <out_dir>/bridge-result.json に書く。
args <- commandArgs(trailingOnly = TRUE)
req <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
out_dir <- normalizePath(req$out_dir, winslash = "/", mustWork = TRUE)
result_path <- file.path(out_dir, "bridge-result.json")

write_result <- function(x) {
  jsonlite::write_json(x, result_path, auto_unbox = TRUE, null = "null", pretty = FALSE)
}
fail <- function(code, message) {
  write_result(list(status = "error", code = code, message = message))
  quit(status = 2, save = "no")
}

if (!requireNamespace("MSLipidMapper", quietly = TRUE)) {
  fail("MSLM_R_PACKAGES_MISSING",
       "R パッケージ MSLipidMapper がありません。上流クローンで R CMD INSTALL を実行してください。")
}
needed <- c("mztab_pathway_measurements", "render_pathway_from_config")
missing_fns <- needed[!needed %in% getNamespaceExports("MSLipidMapper")]
if (length(missing_fns)) {
  fail("MSLM_R_PACKAGES_MISSING",
       paste0("MSLipidMapper が古く、次の関数がありません: ", paste(missing_fns, collapse = ", ")))
}
version <- as.character(utils::packageVersion("MSLipidMapper"))
network_path <- function(p) {
  f <- system.file("extdata", "examples", req$networks[[p]], package = "MSLipidMapper")
  if (!nzchar(f)) fail("PATHWAY_BAD_REQUEST", paste0("同梱パスウェイが見つかりません: ", p))
  f
}

outputs <- list(); mapped <- list()
if (identical(req$mode, "check")) {
  write_result(list(status = "ok", code = NULL, message = "", package_version = version,
                    outputs = outputs, mapped_nodes = mapped))
  quit(status = 0, save = "no")
}

if (identical(req$mode, "measure")) {
  for (p in req$pathways) {
    out_file <- file.path(out_dir, paste0(p, "-measurements.json"))
    res <- tryCatch(
      MSLipidMapper::mztab_pathway_measurements(
        req$mztab, network_path(p), output_file = out_file, group_col = "class",
        acyl_chains = as.character(unlist(req$acyl_chains)), aggregate = req$aggregation,
        normalization = req$normalization, unit = req$unit,
        dataset_id = paste0(req$dataset_id, "-", p)),
      error = function(e) fail(
        if (grepl("acyl", conditionMessage(e), ignore.case = TRUE)) "PATHWAY_BAD_REQUEST" else "MSLM_RUN_FAILED",
        conditionMessage(e)))
    outputs[[p]] <- out_file
    nm <- res$node_mapping
    mapped[[p]] <- if (is.data.frame(nm) && "matched" %in% names(nm)) sum(nm$matched, na.rm = TRUE) else NA
  }
} else if (identical(req$mode, "render")) {
  cfg <- list(
    input = list(format = "mztab-m", file = req$mztab),
    analysis = list(normalization = "none", aggregation = req$aggregation, group_column = "class"),
    pathway = list(browser = req$browser, timeout_seconds = req$timeout_seconds),
    plot = list(type = req$plot_type, x_order = as.list(req$x_order)),
    output = list(paper_size = "A4", orientation = "landscape"))
  cfg_path <- file.path(out_dir, "render-config.yml")
  yaml::write_yaml(cfg, cfg_path)
  for (p in req$pathways) {
    pdf <- file.path(out_dir, paste0(p, ".pdf"))
    res <- tryCatch(
      MSLipidMapper::render_pathway_from_config(cfg_path, network_path = network_path(p),
                                                output_path = pdf),
      error = function(e) fail("PATHWAY_RENDER_FAILED", conditionMessage(e)))
    outputs[[p]] <- pdf
    mapped[[p]] <- unname(unlist(res$mapped_nodes))[1]
  }
} else {
  fail("PATHWAY_BAD_REQUEST", paste0("未知の mode: ", req$mode))
}
write_result(list(status = "ok", code = NULL, message = "", package_version = version,
                  outputs = outputs, mapped_nodes = mapped))
```

`res$node_mapping` の列名（`matched` 等）は上流 `R/pathway_measurements.R` の戻り値を実物で確認して直す（Task 11 Step 3 で実測）。違っていれば `mapped_nodes` は Python 側で JSON-LD から数えた値を正とし、R 側の値は参考値として扱う（Task 8 は JSON-LD 由来の数を出す）。

- [ ] **Step 5: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_bridge.py -q`
Expected: PASS

- [ ] **Step 6: コミット**（バックグラウンド）

```bash
git add metabolomix/pathway/bridge.py metabolomix/pathway/mslm_bridge.R tests/test_pathway_bridge.py
git commit -m "feat(pathway): MSLipidMapper を Rscript で呼ぶブリッジ（タイムアウト・終了判定つき）"
```

---

### Task 7: PDF → PNG

**Files:**
- Create: `metabolomix/plots/pdf_raster.py`
- Test: `tests/test_pdf_raster.py`

**Interfaces:**
- Produces: `MAX_SIDE_PX = 1600`、`pdf_first_page_png(pdf_path: Path, max_side: int = MAX_SIDE_PX) -> tuple[bytes, int, int]`（PNG バイト列, 幅, 高さ）

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pdf_raster.py
import fitz

from metabolomix.plots import pdf_raster


def _pdf(tmp_path, w=842, h=595):          # A4 横（pt）
    doc = fitz.open()
    page = doc.new_page(width=w, height=h)
    page.insert_text((72, 72), "PC")
    path = tmp_path / "x.pdf"
    doc.save(path)
    return path


def test_long_side_is_capped(tmp_path):
    png, w, h = pdf_raster.pdf_first_page_png(_pdf(tmp_path))
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert max(w, h) == 1600
    assert abs(w / h - 842 / 595) < 0.01


def test_small_cap(tmp_path):
    _, w, h = pdf_raster.pdf_first_page_png(_pdf(tmp_path), max_side=400)
    assert max(w, h) == 400
```

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pdf_raster.py -q`
Expected: FAIL（`ImportError`）

- [ ] **Step 3: 実装**

```python
# metabolomix/plots/pdf_raster.py
"""PDF の 1 ページ目を PNG にする（チャット表示用）。

画像トークンは 幅×高さ/750。長辺を上限で揃え、A4 横 1600px ≒ 2,400 トークンに抑える。
"""
from __future__ import annotations

from pathlib import Path

import fitz  # PyMuPDF

MAX_SIDE_PX = 1600


def pdf_first_page_png(pdf_path: Path, max_side: int = MAX_SIDE_PX) -> tuple[bytes, int, int]:
    with fitz.open(pdf_path) as doc:
        page = doc[0]
        rect = page.rect
        zoom = max_side / max(rect.width, rect.height)
        pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        return pix.tobytes("png"), pix.width, pix.height
```

丸めで長辺が 1599 になる場合は、`zoom` を `(max_side + 0.5) / max(...)` にして長辺がちょうど `max_side` になるよう調整する（テストで確認）。

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pdf_raster.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（バックグラウンド）

```bash
git add metabolomix/plots/pdf_raster.py tests/test_pdf_raster.py
git commit -m "feat(plots): PDF の 1 ページ目を長辺上限つきで PNG にする"
```

---

### Task 8: `pathway_project` ツール（セッションスロット・登録）

**Files:**
- Modify: `metabolomix/core/session_state.py`（`PathwayState` と `AnalysisSession.pathway`）
- Create: `metabolomix/pathway/tools.py`
- Modify: `server.py`（`from metabolomix.pathway.tools import *` を `curation_tools` の後に）
- Modify: `tests/test_server_registration.py`（`EXPECTED_TOOLS` に `pathway_project` `pathway_plot`。コメントに経緯を 1 行）
- Test: `tests/test_pathway_tools.py`

**Interfaces:**
- Consumes: Task 1〜7 の全部。既存 `metabolomix.arf.tools._pool_group_labels` / `_require_arf_with_arf2`、`metabolomix.arf2.reader.load_catalog`、`mcp_core._resolve_report_dir`、`mcp_errors.missing_state` / `pathway_error`。
- Produces:
  - `session_state.PathwayState`（`last_run_id: str | None`, `last_run_dir: str | None`）
  - ツール `pathway_project(source="auto", group_a=None, group_b=None, group_a_label="group_a", group_b_label="group_b", pathways=None, acyl_chains=None, aggregation="sum", q_threshold=0.05) -> str`
  - `tools._run_dir(run_id: str) -> Path`、`tools._load_manifest(run_id: str | None) -> dict | str`（`str` は missing_state 封筒。Task 9 が使う）
  - manifest のキー: `run_id source created_at group_a group_b groups pathways acyl_chains aggregation q_threshold normalization dataset_id upstream_commit package_version report files`

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_pathway_tools.py
import json

import numpy as np
import pytest

from metabolomix.core import mcp_core, session_state
from metabolomix.pathway import bridge, tools

from tests.test_pathway_jsonld import _graph  # 小さな JSON-LD 組み立て器を再利用


@pytest.fixture
def lipid_arf(monkeypatch, tmp_path):
    """前処理済み ARF ＋兄弟 .arf2 カタログを持つセッション（ファイル I/O は差し替え）。"""
    s = session_state.AnalysisSession()
    monkeypatch.setattr(session_state, "session", s)
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    s.assay_kind = "lipid"
    a = s.arf
    a.features = [{"MasterAlignmentID": i} for i in (1, 2)]
    a.current_file_path = str(tmp_path / "AlignResult-1.arf")
    a.feature_matrix = np.array([[10.0, 1.0], [11.0, 1.1], [30.0, 0.3], [31.0, 0.31]])
    a.pp_sample_names = ["C_1", "C_2", "P_1", "P_2"]
    a.pp_feature_names = ["Spot_1_height", "Spot_2_height"]
    a.sample_meta = {n: {"role": "sample", "group": n.split("_")[0]} for n in a.pp_sample_names}
    a.preprocessing_recipe = {"normalize": "none", "props": ["height"]}
    catalog = [{"MasterAlignmentID": 1, "Name": "PC 34:1", "Ontology": "PC"},
               {"MasterAlignmentID": 2, "Name": "LPC 18:1", "Ontology": "LPC"}]
    monkeypatch.setattr(tools, "_arf2_catalog", lambda: (catalog, None))
    return s


@pytest.fixture
def fake_measure(monkeypatch):
    """run_bridge を差し替え、measure で JSON-LD を書く。"""
    calls = []

    def run(request, run_dir, *, timeout_s=300.0):
        calls.append(request)
        out = {}
        for p in request["pathways"]:
            path = run_dir / f"{p}-measurements.json"
            doc = _graph()
            # 群名を実際の比較に合わせる
            for x in doc["@graph"]:
                if x["@type"] == "mslm:Sample":
                    x["sampleGroup"] = "C" if x["identifier"] == "S1" else "P"
            path.write_text(json.dumps(doc), encoding="utf-8")
            out[p] = str(path)
        return {"status": "ok", "outputs": out, "mapped_nodes": {}, "package_version": "0.1"}

    monkeypatch.setattr(bridge, "run_bridge", run)
    monkeypatch.setattr(bridge, "upstream_commit", lambda: "abc123")
    return calls


def test_requires_lipid_assay(lipid_arf, fake_measure):
    lipid_arf.assay_kind = "metabolite"
    body = json.loads(tools.pathway_project(group_a="C", group_b="P"))
    assert body["error"]["code"] == "missing_state"
    assert body["error"]["required_tools"] == ["record_objective"]


def test_requires_preprocessed_matrix(lipid_arf, fake_measure):
    lipid_arf.arf.feature_matrix = None
    body = json.loads(tools.pathway_project(group_a="C", group_b="P"))
    assert body["error"]["required_tools"] == ["arf_preprocess"]


def test_unknown_pathway_is_bad_request(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["nope"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_project_writes_run_and_returns_tables(lipid_arf, fake_measure, tmp_path):
    body = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))
    run_id = body["run_id"]
    run_dir = tmp_path / "reports" / "pathway" / run_id
    assert (run_dir / "input.mztab").is_file()
    assert (run_dir / "nodes.tsv").is_file()
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["upstream_commit"] == "abc123"
    assert manifest["groups"] == {"C": 2, "P": 2}
    assert body["pathways"]["remodeling"]["nodes_tsv"].startswith(
        "node\tscope\tn_features\tmean_a\tmean_b\tlog2fc\tq\tcall")
    assert body["pathways"]["remodeling"]["unmatched_subclasses"] == ["TG"]
    assert "not curated reactions" in body["notes"]["adjacency"] or "反応" in body["notes"]["adjacency"]
    assert lipid_arf.pathway.last_run_id == run_id
    req = fake_measure[0]
    assert req["mode"] == "measure" and req["networks"]["remodeling"] == "remodeling_lipidonly.cyjs"


def test_bridge_error_becomes_envelope(lipid_arf, monkeypatch):
    def boom(*a, **k):
        raise bridge.BridgeError("MSLM_NOT_CONFIGURED", "設定して", {"setting": "x"})
    monkeypatch.setattr(bridge, "run_bridge", boom)
    body = json.loads(tools.pathway_project(group_a="C", group_b="P"))
    assert body["error"]["code"] == "MSLM_NOT_CONFIGURED"
```

`_graph()` の JSON-LD は S1/S2 の 2 サンプルしか持たないので検定値は NaN（`call=NA`）になるが、このテストは表の形と保存物を見るだけでよい。数値の正しさは Task 5 が持つ。

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py -q`
Expected: FAIL（`ImportError: metabolomix.pathway.tools`）

- [ ] **Step 3: セッションスロットを足す**

`metabolomix/core/session_state.py` の `CurationState` の後に:

```python
class PathwayState:
    """パスウェイ投影（pathway_project / pathway_plot）の直近 run へのポインタ。

    run の正準は `reports/pathway/<run_id>/manifest.json`。ここは run_id 省略時の
    既定を決めるためだけにある。他パーサのスロットに触らない。
    """

    def __init__(self):
        self.last_run_id: str | None = None
        self.last_run_dir: str | None = None
```

`AnalysisSession.__init__` の `self.curation = CurationState()` の後に:

```python
        # --- パスウェイ投影スロット（MSLipidMapper。他スロットとは独立） ---
        self.pathway = PathwayState()
```

- [ ] **Step 4: ツールを実装**

```python
# metabolomix/pathway/tools.py
"""MCP ツール: MSLipidMapper を通したパスウェイ投影（spec 2026-10-14）。"""
from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core import mcp_core, mcp_errors, session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.pathway import bridge, compare, handoff, jsonld

__all__ = ["pathway_project"]

_NODE_COLUMNS = ["node", "scope", "n_features", "mean_a", "mean_b", "log2fc", "q", "call"]
_PAIR_COLUMNS = ["node_1", "call_1", "log2fc_1", "node_2", "call_2", "log2fc_2", "via"]
_ADJACENCY_NOTE = (
    "隣接は MSLipidMapper のパスウェイ図のつながり（接合ノードだけを挟むもの）で、"
    "上流は反応だとは主張していない（not curated reactions）。逆向きの組は観察の提示であり、"
    "基質・生成物の解釈は文献で裏を取ること。"
    "掲載条件は両端とも q < 閾値かつ向きが逆。片方が ns/NA の組は件数だけ示す。")


def _arf2_catalog():
    """兄弟 .arf2 のスポット一覧。無ければ (None, missing_state 封筒)。"""
    from metabolomix.arf.tools import _require_arf_with_arf2
    from metabolomix.arf2.reader import load_catalog
    found = _require_arf_with_arf2()
    if isinstance(found, str):
        return None, found
    return load_catalog(found), None


def _run_dir(run_id: str) -> Path:
    return mcp_core._resolve_report_dir() / "pathway" / run_id


def _load_manifest(run_id: str | None):
    run_id = run_id or session_state.session.pathway.last_run_id
    if not run_id:
        return mcp_errors.missing_state("pathway_run", ["pathway_project"],
                                        "先に pathway_project を実行してください。")
    path = _run_dir(run_id) / "manifest.json"
    if not path.is_file():
        return mcp_errors.missing_state("pathway_run", ["pathway_project"],
                                        f"run {run_id} の manifest.json がありません。")
    return json.loads(path.read_text(encoding="utf-8"))


def _pick_source(source: str) -> str:
    s = session_state.session
    if source != "auto":
        return source
    if s.arf.feature_matrix is not None:
        return "arf"
    if s.dataset is not None:
        return "dataset"
    return "arf"


def _recipe_text(recipe: dict) -> str:
    parts = [f"{k}={v}" for k, v in sorted((recipe or {}).items()) if v not in (None, "", [])]
    return "; ".join(parts) or "none"


def _build_arf(group_a, group_b):
    from metabolomix.arf.tools import _pool_group_labels
    a = session_state.session.arf
    if a.feature_matrix is None:
        return None, mcp_errors.missing_state(
            "preprocessed_matrix", ["arf_preprocess"],
            "先に arf_preprocess を実行してください（前処理後行列が必要）。")
    catalog, err = _arf2_catalog()
    if err:
        return None, err
    labels = [None if (a.sample_meta.get(n) or {}).get("role", "sample") != "sample"
              else (a.sample_meta.get(n) or {}).get("group") for n in a.pp_sample_names]
    try:
        _, _, resolved = _pool_group_labels(a.pp_sample_names, labels, group_a, group_b)
    except ValueError as exc:
        raise handoff.HandoffError("PATHWAY_BAD_REQUEST", str(exc)) from exc
    table = handoff.build_from_arf(
        matrix=a.feature_matrix, sample_names=a.pp_sample_names,
        feature_names=a.pp_feature_names, sample_meta=a.sample_meta,
        excluded_samples=a.excluded_samples, excluded_spots=a.excluded_spots,
        catalog={int(s["MasterAlignmentID"]): s for s in catalog},
        props=(a.preprocessing_recipe or {}).get("props"),
        group_a=str(group_a), group_b=str(group_b), compared=resolved)
    return (table, _recipe_text(a.preprocessing_recipe), a.dataset_id,
            str(group_a), str(group_b)), None


def _build_dataset(group_a, group_b, label_a, label_b):
    ds = session_state.session.dataset
    if ds is None:
        return None, mcp_errors.missing_state("dataset", ["dataset_load"],
                                              "先に dataset_load を実行してください。")
    if ds.pp_matrix is None:
        return None, mcp_errors.missing_state("preprocessed_matrix", ["dataset_preprocess"],
                                              "先に dataset_preprocess を実行してください。")
    table = handoff.build_from_dataset(
        pp_matrix=ds.pp_matrix, pp_sample_names=ds.pp_sample_names,
        pp_feature_names=ds.pp_feature_names, roles=ds.roles,
        sample_names=ds.sample_names, sample_ids=ds.sample_ids,
        metadata_rows=ds.sample_metadata_rows, annotations=ds.feature_annotations,
        group_a=list(group_a or []), group_b=list(group_b or []),
        group_a_label=label_a, group_b_label=label_b)
    return (table, _recipe_text(ds.preprocessing_recipe), ds.dataset_id,
            label_a, label_b), None


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, openWorldHint=False),
          structured_output=False)
def pathway_project(source: str = "auto", group_a=None, group_b=None,
                    group_a_label: str = "group_a", group_b_label: str = "group_b",
                    pathways: list[str] | None = None, acyl_chains: list[str] | None = None,
                    aggregation: str = "sum", q_threshold: float = 0.05) -> str:
    """前処理済みの脂質データを MSLipidMapper のパスウェイへ投影し、2 群を比べる。

    返すのは小さな表だけ（ノード表と、逆向きに動く隣接ペア）。JSON-LD・中間 mzTab-M は
    解析フォルダの reports/pathway/<run_id>/ に残る。図は pathway_plot。

    - source: auto | arf | dataset。arf は arf_preprocess 後、dataset は dataset_preprocess 後。
    - group_a / group_b: arf は Class ID か因子トークン（arf_differential と同じ）、
      dataset はサンプル名のリスト（group_a_label / group_b_label が群名になる）。
    - pathways: remodeling / ceramide / global（既定は 3 つ全部）。
    - acyl_chains: 例 ["18:0","20:4"]。指定すると鎖サブセットのスコープも比べる。
    - aggregation: sum | mean | median（サブクラス内の分子種の集計）。
    検定は log2(x+1) の Welch、BH はパスウェイ×スコープ種別ごと。log2fc は正なら group_b が高い。
    脂質アッセイ（record_objective の assay_kind=lipid）でだけ使える。
    """
    s = session_state.session
    if s.assay_kind != "lipid":
        return mcp_errors.missing_state(
            "assay_kind", ["record_objective"],
            "MSLipidMapper は脂質専用です。record_objective(assay_kind='lipid') で種別を確定してください。")
    pathways = list(pathways or bridge.PATHWAY_NETWORKS)
    unknown = [p for p in pathways if p not in bridge.PATHWAY_NETWORKS]
    if unknown:
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST", f"未知のパスウェイ: {', '.join(unknown)}",
            {"valid": list(bridge.PATHWAY_NETWORKS)})
    if aggregation not in ("sum", "mean", "median"):
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST",
                                        "aggregation は sum / mean / median のどれかです。")
    if group_a is None or group_b is None:
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST", "group_a と group_b を指定してください。")

    src = _pick_source(source)
    try:
        built, err = (_build_arf(group_a, group_b) if src == "arf"
                      else _build_dataset(group_a, group_b, group_a_label, group_b_label))
    except handoff.HandoffError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message, exc.details)
    if err:
        return err
    table, normalization, dataset_id, label_a, label_b = built

    run_id = f"pw_{datetime.now():%Y%m%d_%H%M%S}_{uuid.uuid4().hex[:6]}"
    run_dir = _run_dir(run_id)
    run_dir.mkdir(parents=True, exist_ok=True)
    mztab_path = run_dir / "input.mztab"
    handoff.write_mztab(table, mztab_path, normalization=normalization)
    request = {"mode": "measure", "mztab": str(mztab_path), "out_dir": str(run_dir),
               "pathways": pathways,
               "networks": {p: bridge.PATHWAY_NETWORKS[p] for p in pathways},
               "acyl_chains": list(acyl_chains or []), "aggregation": aggregation,
               "normalization": normalization, "unit": "arbitrary", "dataset_id": run_id}
    try:
        result = bridge.run_bridge(request, run_dir)
        graphs = {p: jsonld.load_graph(Path(result["outputs"][p])) for p in pathways}
    except bridge.BridgeError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message, exc.details)
    except jsonld.JsonldError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message)

    per_pathway, all_rows = {}, []
    for p, graph in graphs.items():
        rows = compare.compare_scopes(graph, label_a, label_b, q_threshold)
        pairs, undeterminable = compare.opposite_pairs(graph, rows)
        lipid_nodes = sum(1 for k in graph.node_kinds.values() if k != "junction")
        valued = len({r["node_id"] for r in rows if r["scope"] == "whole"})
        per_pathway[p] = {
            "mapped_nodes": f"{valued}/{lipid_nodes}",
            "unmatched_subclasses": graph.unmatched_subclasses,
            "nodes_tsv": compare.rows_to_tsv(rows, _NODE_COLUMNS),
            "opposite_pairs_tsv": compare.rows_to_tsv(pairs, _PAIR_COLUMNS),
            "undeterminable_pairs": undeterminable,
        }
        all_rows += [dict(r, pathway=p) for r in rows]
    (run_dir / "nodes.tsv").write_text(
        compare.rows_to_tsv(all_rows, ["pathway", *_NODE_COLUMNS]), encoding="utf-8")

    groups = {}
    for g in table.groups:
        groups[g] = groups.get(g, 0) + 1
    manifest = {"run_id": run_id, "source": src,
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "group_a": label_a, "group_b": label_b, "groups": groups,
                "pathways": pathways, "acyl_chains": list(acyl_chains or []),
                "aggregation": aggregation, "q_threshold": q_threshold,
                "normalization": normalization, "dataset_id": dataset_id,
                "upstream_commit": bridge.upstream_commit(),
                "package_version": result.get("package_version"),
                "report": table.report,
                "files": {"mztab": "input.mztab", "nodes": "nodes.tsv",
                          **{p: Path(result["outputs"][p]).name for p in pathways}}}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1),
                                           encoding="utf-8")
    s.pathway.last_run_id = run_id
    s.pathway.last_run_dir = str(run_dir)

    return json_payload({
        "run_id": run_id, "source": src, "run_dir": str(run_dir),
        "groups": groups, "group_a": label_a, "group_b": label_b,
        "normalization": normalization,
        "handoff": table.report,
        "pathways": per_pathway,
        "notes": {"adjacency": _ADJACENCY_NOTE,
                  "values": "MSLipidMapper は正規化しない。値はこちらの前処理後行列をサブクラスごとに集計したもの。",
                  "figure": "図は pathway_plot(run_id) で表示できる。"},
    })
```

テストの `"not curated reactions" in ...` は `_ADJACENCY_NOTE` に `not curated reactions` が含まれることで満たす（上の文面のとおり）。`ToolAnnotations` は既存の書き込み系ツール（`write_report` 等）と同じ値を使う——実装前に `grep -n "def write_report" -B2 metabolomix/tools/reports.py` で合わせる。

`server.py` に追加（`curation_tools` の行の後）:

```python
from metabolomix.pathway.tools import *  # pathway_project, pathway_plot（MSLipidMapper 連携）
```

`tests/test_server_registration.py` の `EXPECTED_TOOLS` に `"pathway_project"` を足し、冒頭の経緯コメントに `# パスウェイ投影: pathway_project / pathway_plot を追加(74→76)。` を追記（`pathway_plot` は Task 9 で足す）。

- [ ] **Step 5: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py tests/test_server_registration.py -q`
Expected: PASS（`test_readme_links` と `test_workflow_docs` はこの時点で落ちる。Task 10 で直す）

- [ ] **Step 6: コミット**（バックグラウンド。フックが全テストを走らせるので、Task 10 まで `--no-verify` を使わず、ここでは Task 9・10 と合わせて 1 コミットにしてもよい。その場合このステップを飛ばし、Task 10 の最後でまとめてコミットする）

---

### Task 9: `pathway_plot` ツール

**Files:**
- Modify: `metabolomix/pathway/tools.py`
- Modify: `tests/test_server_registration.py`（`"pathway_plot"`）
- Test: `tests/test_pathway_tools.py`（追記）

**Interfaces:**
- Consumes: Task 8 の `_load_manifest` / `_run_dir`、Task 6 の `bridge.run_bridge` / `browser_path`、Task 7 の `pdf_raster.pdf_first_page_png`。
- Produces: ツール `pathway_plot(run_id=None, pathways=None, plot_type="dot") -> list`（先頭が説明文 str、続いて `Image`）または封筒 str。

- [ ] **Step 1: 失敗するテストを書く**（`tests/test_pathway_tools.py` に追記）

```python
import fitz
from mcp.server.fastmcp import Image


def _write_pdf(path):
    doc = fitz.open()
    doc.new_page(width=842, height=595).insert_text((72, 72), "PC")
    doc.save(path)


def test_plot_without_run_is_missing_state(lipid_arf):
    body = json.loads(tools.pathway_plot())
    assert body["error"]["required_tools"] == ["pathway_project"]


def test_plot_renders_pngs(lipid_arf, fake_measure, monkeypatch):
    run_id = json.loads(tools.pathway_project(group_a="C", group_b="P",
                                              pathways=["remodeling", "ceramide"]))["run_id"]
    seen = []

    def render(request, run_dir, *, timeout_s=300.0):
        seen.append(request)
        out = {}
        for p in request["pathways"]:
            pdf = run_dir / f"{p}.pdf"
            _write_pdf(pdf)
            out[p] = str(pdf)
        return {"status": "ok", "outputs": out, "mapped_nodes": {}}

    monkeypatch.setattr(bridge, "run_bridge", render)
    out = tools.pathway_plot(pathways=["ceramide"])
    assert isinstance(out, list) and isinstance(out[0], str) and isinstance(out[1], Image)
    assert len(out) == 2
    assert "ceramide.pdf" in out[0] and "q" in out[0]
    assert seen[0]["mode"] == "render" and seen[0]["x_order"] == ["C", "P"]
    assert seen[0]["mztab"].endswith("input.mztab")
    assert run_id in seen[0]["out_dir"]


def test_plot_pathway_not_in_run_is_bad_request(lipid_arf, fake_measure):
    json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))
    body = json.loads(tools.pathway_plot(pathways=["global"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_render_failure_envelope(lipid_arf, fake_measure, monkeypatch):
    json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))

    def fail(*a, **k):
        raise bridge.BridgeError("PATHWAY_RENDER_FAILED", "Chromium がありません")
    monkeypatch.setattr(bridge, "run_bridge", fail)
    body = json.loads(tools.pathway_plot())
    assert body["error"]["code"] == "PATHWAY_RENDER_FAILED"
```

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py -q -k plot`
Expected: FAIL（`AttributeError: pathway_plot`）

- [ ] **Step 3: 実装**（`metabolomix/pathway/tools.py` に追記し、`__all__` を `["pathway_project", "pathway_plot"]` に）

```python
from mcp.server.fastmcp import Image

from metabolomix.plots import pdf_raster

_FIGURE_NOTE = ("図のノード内の p 値ラベルは MSLipidMapper（Shiny と同じ計算）のもので、"
                "pathway_project の q 値とは別物。")


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, openWorldHint=True),
          structured_output=False)
def pathway_plot(run_id: str | None = None, pathways: list[str] | None = None,
                 plot_type: str = "dot"):
    """pathway_project の run をパスウェイ図（PDF）に描き、チャットに PNG で表示する。

    - run_id: 省略時は直近の run。
    - pathways: その run に含まれるものから選ぶ（既定は全部）。1 枚約 2,400 画像トークン。
    - plot_type: dot | box | violin（ノード内のグラフ）。
    PDF 描画にはヘッドレス Chromium（[mslipidmapper] browser）と unpkg.com への接続が要る。
    """
    manifest = _load_manifest(run_id)
    if isinstance(manifest, str):
        return manifest
    chosen = list(pathways or manifest["pathways"])
    missing = [p for p in chosen if p not in manifest["pathways"]]
    if missing:
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST", f"この run に無いパスウェイ: {', '.join(missing)}",
            {"run_pathways": manifest["pathways"]})
    if plot_type not in ("dot", "box", "violin"):
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST",
                                        "plot_type は dot / box / violin のどれかです。")
    run_dir = _run_dir(manifest["run_id"])
    others = [g for g in manifest["groups"] if g not in (manifest["group_a"], manifest["group_b"])]
    request = {"mode": "render", "mztab": str(run_dir / manifest["files"]["mztab"]),
               "out_dir": str(run_dir), "pathways": chosen,
               "networks": {p: bridge.PATHWAY_NETWORKS[p] for p in chosen},
               "aggregation": manifest["aggregation"], "browser": bridge.browser_path(),
               "plot_type": plot_type,
               "x_order": [manifest["group_a"], manifest["group_b"], *others],
               "timeout_seconds": 90}
    try:
        result = bridge.run_bridge(request, run_dir, timeout_s=120.0 * len(chosen))
    except bridge.BridgeError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message, exc.details)

    images, lines = [], [f"run {manifest['run_id']}（{manifest['group_a']} 対 {manifest['group_b']}）"]
    for p in chosen:
        pdf = Path(result["outputs"][p])
        png, w, h = pdf_raster.pdf_first_page_png(pdf)
        images.append(Image(data=png, format="png"))
        lines.append(f"- {p}: {pdf}（{w}×{h}px）")
    lines.append(_FIGURE_NOTE)
    return ["\n".join(lines), *images]
```

`tests/test_server_registration.py` の `EXPECTED_TOOLS` に `"pathway_plot"` を足す。

- [ ] **Step 4: 通ることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py tests/test_server_registration.py -q`
Expected: PASS

---

### Task 10: 文書（USAGE・workflow・output_format・CLAUDE.md）

**Files:**
- Modify: `USAGE.md`（見出し件数 `全74ツール` → `全76ツール`、新節 `## 17. パスウェイ投影(MSLipidMapper)`）
- Create: `docs/workflow/pathway.md`
- Modify: `docs/workflow/index.md`（目次に行を足し、合計を 55 / 総数 76 に）
- Modify: `tests/test_workflow_docs.py`（`IN_SCOPE["pathway.md"] = ("pathway_project", "pathway_plot")`）
- Create: `docs/output_format/pathway.md`
- Modify: `metabolomix/core/mcp_core.py`（`OUTPUT_FORMAT_SECTIONS["pathway"]`）
- Modify: `metabolomix/tools/resources.py`（2 つのリソース description のトピック列挙に `pathway`）
- Modify: `docs/output_format/core.md`（トピック索引に `pathway` があればそこへ 1 行）
- Modify: `CLAUDE.md`（冒頭「ツール 74」→「ツール 76」、構成図に `metabolomix/pathway/` の 1 行、読む形式に JSON-LD の 1 文）

- [ ] **Step 1: 文書テストが落ちていることを確認**

Run: `C:/Python314/python.exe -m pytest tests/test_readme_links.py tests/test_workflow_docs.py -q`
Expected: FAIL（USAGE 件数・CLAUDE.md 規模表記・workflow の対象範囲）

- [ ] **Step 2: USAGE.md**

見出しを `# USAGE — ms-data-parser MCP ツール一覧(全76ツール)` にし、`## 16. サーバ自身の更新` の前に挿入（以降の節番号は既存の付け方に合わせて振り直す）:

```markdown
## 16. パスウェイ投影(MSLipidMapper)

前処理済みの脂質データを MSLipidMapper(ローカルの Rscript)で remodeling / ceramide / global の
パスウェイへ投影する。`lipidmix.local.toml` の `[mslipidmapper] rscript`(と描画用の `browser`)が要る。
脂質アッセイ(`record_objective(assay_kind="lipid")`)でだけ使える。

| ツール | 機能 |
|--------|------|
| `pathway_project` | `arf_preprocess` / `dataset_preprocess` 後の行列から最小 mzTab-M を書き、MSLipidMapper でサブクラス(と任意のアシル鎖サブセット)をノードへ集計する。2 群を log2 空間の Welch＋BH(パスウェイ×スコープ種別ごと)で比べ、ノード表と「逆向きに動く隣接ペア」を TSV で返す。JSON-LD・mzTab・manifest は解析フォルダの `reports/pathway/<run_id>/`。 |
| `pathway_plot` | `pathway_project` の run を MSLipidMapper で PDF に描き、1 ページ目を PNG(長辺 1600px)でチャットに表示する。ヘッドレス Chromium と unpkg.com への接続が要る。図中の p 値は MSLipidMapper のもの。 |
```

- [ ] **Step 3: docs/workflow/pathway.md**

```markdown
# ワークフロー: パスウェイ投影（MSLipidMapper）

前処理済みの脂質データを、ローカルの Rscript で MSLipidMapper（R パッケージ）に渡し、
同梱パスウェイ（remodeling / ceramide / global）へ投影する経路。前提状態の連鎖は
`record_objective(assay_kind="lipid")` ＋ `arf_preprocess`（または `dataset_preprocess`）→
`pathway_project` → `pathway_plot`。値の**意味**は
[pathway.md（output_format）](../output_format/pathway.md) を見ること。

`session.pathway`（他スロットとは独立）は直近の `run_id` だけを持つ。run の正準は
`reports/pathway/<run_id>/manifest.json`。

## pathway_project

前提: `session.assay_kind == "lipid"`（違えば `missing_state` → `record_objective`）。
arf は `session.arf.feature_matrix` と兄弟 `.arf2`、dataset は `DatasetState.pp_matrix`。
状態変更: `session.pathway.last_run_id` / `last_run_dir`。

1. metabolomix/pathway/tools.py  pathway_project()
2. ├─ metabolomix/pathway/tools.py  _build_arf()
3. │  ├─ metabolomix/arf/tools.py  _pool_group_labels()
4. │  ├─ metabolomix/pathway/tools.py  _arf2_catalog()
5. │  └─ metabolomix/pathway/handoff.py  build_from_arf()
6. ├─ metabolomix/pathway/tools.py  _build_dataset()
7. │  └─ metabolomix/pathway/handoff.py  build_from_dataset()
8. ├─ metabolomix/pathway/handoff.py  write_mztab()
9. ├─ metabolomix/pathway/bridge.py  run_bridge()
10. │  └─ metabolomix/pathway/mslm_bridge.R（mode=measure → MSLipidMapper::mztab_pathway_measurements）
11. ├─ metabolomix/pathway/jsonld.py  load_graph()
12. ├─ metabolomix/pathway/compare.py  compare_scopes()
13. │  └─ metabolomix/analysis/differential.py  two_group_test()
14. └─ metabolomix/pathway/compare.py  opposite_pairs()
15.    └─ metabolomix/pathway/jsonld.py  adjacent_lipid_pairs()

## pathway_plot

前提: `pathway_project` の run（`run_id` 省略時は `session.pathway.last_run_id`）。
無ければ `missing_state` → `pathway_project`。状態変更: なし（PDF を run ディレクトリに書く）。

1. metabolomix/pathway/tools.py  pathway_plot()
2. ├─ metabolomix/pathway/tools.py  _load_manifest()
3. ├─ metabolomix/pathway/bridge.py  run_bridge()
4. │  └─ metabolomix/pathway/mslm_bridge.R（mode=render → MSLipidMapper::render_pathway_from_config）
5. └─ metabolomix/plots/pdf_raster.py  pdf_first_page_png()
```

`tests/test_workflow_docs.py` の参照抽出が `.R` ファイル名や括弧書きを誤検出する場合は、
10・4 行目を「（R 側。`mslm_bridge.R` が上流関数を呼ぶ）」のような地の文に変えてテストを通す
（既存文書の書式に合わせる。実装前に `tests/test_workflow_docs.py` の参照正規表現を読む）。

- [ ] **Step 4: docs/workflow/index.md**

目次表の `curation.md` の行の後に
`| [pathway.md](pathway.md) | パスウェイ投影（MSLipidMapper）— サブクラスのノード集計・2 群比較・図 | 2 |`
を足し、`合計 53 ツール` → `合計 55 ツール`、`登録ツール総数は 74（53 + 21）` → `登録ツール総数は 76（55 + 21）`。

- [ ] **Step 5: docs/output_format/pathway.md**

```markdown
# 出力形式: パスウェイ投影（pathway_project / pathway_plot）

## 何の値か

- MSLipidMapper は**正規化しない**。値はこちらの前処理後行列（`normalization` 欄にレシピ）を、
  パスウェイノードに対応するサブクラスごとに `aggregation`（既定 sum）で集計したもの。
- 中間 mzTab-M に載るのは生体試料（role=sample）だけ。サブクラスが空・Unknown・Others の行と、
  手動除外スポットは載らない（件数は `handoff`）。同名の行は平均強度最大の 1 行だけ残す
  （`dropped_duplicate_names`）——sum がアダクト違いを二重計上しないため。
- arf はレシピ先頭の prop（`handoff.prop`）だけを使う。

## ノード表（nodes_tsv）

| 列 | 意味 |
|---|---|
| node | パスウェイのノード名（＝MS-DIAL Ontology と完全一致したサブクラス） |
| scope | `whole`（サブクラス全体）か、アシル鎖（`18:0` 等。そのノードの鎖サブセット） |
| n_features | そのスコープに入った分子種の数 |
| mean_a / mean_b | 各群の集計値の平均（生スケール） |
| log2fc | log2(x+1) 空間の群平均差。**正なら group_b が高い** |
| q | BH-FDR。補正単位は**パスウェイ × スコープ種別**（whole と各鎖は別々） |
| call | up / down（q < 閾値）、ns、NA（片群の値が足りない等で検定不能） |

## 逆向きに動く隣接ペア（opposite_pairs_tsv）

- **隣接**＝間に接合ノード（PathwayJunction）だけを挟んでつながる 2 つの脂質ノード（`via` が
  direct / junction）。値の付かない化合物ノード（CDP-DAG 等）を挟むものは隣接としない。
- 掲載は**両端とも q < 閾値かつ向きが逆**の組だけ。片方が ns/NA の組は `undeterminable_pairs` に件数のみ。
- MSLipidMapper のつながりは描画用テンプレートで、反応の主張ではない（上流文書:
  "not automatically curated reactions"）。**観察の提示**として扱い、基質・生成物の解釈は文献で裏を取る。

## 対応付けの見出し

- `mapped_nodes` = 値の付いたノード / 接合ノード以外の全ノード。
- `unmatched_subclasses` = データにあるがどのノードにも乗らなかったサブクラス
  （ノード名との完全一致のみ。`EtherPC` などの表記揺れはここに出る）。

## 図（pathway_plot）

- 各パスウェイの PDF（`reports/pathway/<run_id>/<pathway>.pdf`）の 1 ページ目を PNG にしたもの。
- ノード内の p 値ラベルは MSLipidMapper（Shiny と同じ計算）のもので、ノード表の q とは**別物**。
```

`metabolomix/core/mcp_core.py` の `OUTPUT_FORMAT_SECTIONS` に
`"pathway": "パスウェイ投影（pathway_project / pathway_plot、MSLipidMapper 連携）",` を足す。
`metabolomix/tools/resources.py` の 2 つの description の列挙に `pathway` を足す
（`"... library, curation."` → `"... library, curation, pathway."`、テンプレート側は
`"curation (...), pathway (MSLipidMapper pathway projection: node table and adjacent opposite pairs), or core ..."`）。
`MCP_INSTRUCTIONS` の topic 列挙（`library`, or `curation` …）にも `pathway` を同じ形で足す
（`grep -n "or \`curation\`" metabolomix/core/mcp_core.py`）。

- [ ] **Step 6: CLAUDE.md**

- 冒頭「ツール 74・リソース 5・リソーステンプレート 3」→「ツール 76・…」。
- 「読む形式」の段落末尾に 1 文: 「MSLipidMapper（R）が書く**パスウェイ測定 JSON-LD**も読む（`pathway_project`）。」
- 構成図の `metabolomix/tools/` の行の前に:
  `metabolomix/pathway/   MSLipidMapper 連携（最小 mzTab-M の受け渡し・Rscript ブリッジ・JSON-LD 読み取り・隣接比較）`
- 「外部資産の場所」の段落の設定キー列挙に `[mslipidmapper] rscript` `browser` `repo`、環境変数 `MSLIPIDMAPPER_RSCRIPT` `MSLIPIDMAPPER_BROWSER` `MSLIPIDMAPPER_REPO` を足す。

- [ ] **Step 7: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全部 PASS（R・Rscript 無しで）

- [ ] **Step 8: コミット**（Task 8・9 をまとめる場合はここで。バックグラウンド）

```bash
git add metabolomix/core/session_state.py metabolomix/pathway/tools.py server.py \
  tests/test_server_registration.py tests/test_pathway_tools.py USAGE.md \
  docs/workflow/pathway.md docs/workflow/index.md tests/test_workflow_docs.py \
  docs/output_format/pathway.md metabolomix/core/mcp_core.py metabolomix/tools/resources.py \
  docs/output_format/core.md CLAUDE.md
git commit -m "feat(pathway): pathway_project / pathway_plot ツールと文書"
```

---

### Task 11: R 環境の準備と実データ検証（手動・ユーザー同席）

**Files:**
- Modify（追跡外）: `lipidmix.local.toml`、`docs/HISTRY.md`、`docs/task.md`
- Modify（リポジトリ外）: `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`
- Modify（必要なら）: `docs/superpowers/specs/2026-09-04-mslipidmapper-transferable-assets.md` §4-2 に実測結果を追記

- [ ] **Step 1: R 依存と MSLipidMapper をインストールする**（**ユーザーに確認してから**。CRAN/Bioconductor から大量に入る）

```bash
"/c/Program Files/R/R-4.6.0/bin/Rscript" -e 'install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("SummarizedExperiment","S4Vectors","clusterProfiler","ComplexHeatmap"), ask=FALSE, update=FALSE)'
"/c/Program Files/R/R-4.6.0/bin/Rscript" -e 'd <- read.dcf("C:/Users/yuu18/MSLipidMapper/DESCRIPTION", "Imports"); p <- trimws(strsplit(gsub("\\(.*?\\)", "", d), ",")[[1]]); p <- setdiff(p[nzchar(p)], rownames(installed.packages())); if (length(p)) BiocManager::install(p, ask=FALSE, update=FALSE)'
"/c/Program Files/R/R-4.6.0/bin/R" CMD INSTALL --no-test-load "C:/Users/yuu18/MSLipidMapper"
```

ユーザーライブラリへの書き込みになる。失敗したパッケージ名を記録し、`docs/HISTRY.md` に残す。

- [ ] **Step 2: 設定ファイルを書く**（`lipidmix.local.toml`。追跡外）

```toml
[mslipidmapper]
rscript = 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe'
browser = 'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
repo = 'C:\Users\yuu18\MSLipidMapper'
```

- [ ] **Step 3: check モードとノード対応の実測**

`check.py`（スクラッチ）で `bridge.run_bridge({"mode": "check", "out_dir": <tmp>}, <tmp>)` を呼び、`package_version` が返ることを確かめる。続けて `mztab_pathway_measurements` の戻り値 `node_mapping` の列名を R で確認し（`str(res$node_mapping)`）、`mslm_bridge.R` の `matched` 参照を実物に合わせる（違っていれば修正をコミット）。

- [ ] **Step 4: 実データで pathway_project / pathway_plot**（MCP 経由。neg と pos それぞれ）

`load_dataset("C:\Users\yuu18\datasets\20260930_EV\neg_wiff")` → `record_objective(assay_kind="lipid", ...)`
→ `arf_preprocess(...)` → `pathway_project(group_a="cell_C", group_b="cell_PlnA")`
→ `pathway_plot()`。pos も同様。記録すること:
  - 各パスウェイの `mapped_nodes` と `unmatched_subclasses`（spec 2026-09-04 §4-2 の答え）
  - 逆向き隣接ペアの件数、判定不能件数
  - 所要時間（measure / render）
  - PNG がチャットに出たか（Claude Code / Claude Desktop）
  群名は `arf_list_classes` の実物に合わせる。

- [ ] **Step 5: 記録**

- `docs/HISTRY.md` に日付見出しで追記（実測値・インストールで詰まった点）。
- `docs/task.md` に残課題があれば追記（例: `unmatched_subclasses` の語彙対応表）。
- vault の流れ図に「脂質アッセイ → 前処理 → `pathway_project` → `pathway_plot`」の枝を足し、末尾の出典行の日付を更新。
- `C:/Users/yuu18/.claude/kb` に再利用できる知見（R の依存インストール手順・上流の JSON-LD の形）があれば INDEX.md を見て重複を確認のうえ 1 エントリ。

- [ ] **Step 6: 仕上げ**

`superpowers:finishing-a-development-branch` に従い、`main` へ `--no-ff` マージ（`Merge feat/mslm-pathway: MSLipidMapper を通したパスウェイ投影ツールを足す`）。push はユーザーに確認。
