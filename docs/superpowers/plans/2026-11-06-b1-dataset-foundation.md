# B-1 データセットの土台（解析フォルダのアダプタ・照合・ハンドルの一本化・ライブラリの記憶）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata` / `_tags.xml` / `.dcl` / `.EIC.aef`）から `DatasetState` を作るアダプタを足し、mzTab-M のアダプタとの照合で同じデータになることを縛り、読み込んだデータセットをすべてハンドルの登録簿に入れ、データセットが自分を同定したライブラリを覚えるようにする。

**Architecture:** 形式に依存しないデータセットの組み立てを新しいパッケージ `metabolomix/dataset/` に置く（ファイルの選択 `alignment_files`・アダプタ `msdial_adapter` と同定の写し `msdial_identity`・ライブラリの記憶 `library_memory`・照合 `agreement`・入口 `opening`）。`.arf2` の欄の読み取りは `metabolomix/arf2/spot_fields.py`、`.mddata` は `metabolomix/msdial/classes.py` に足す。公開面は**今のツール名のまま**で `dataset_load` と `dataset_handle_open` に `directory` / `source` / `measure` を足すだけにする（改名・旧ツールの削除は B-2）。対話のセッションの `dataset` / `library` は「直近のハンドル」を指す property にする。

**Tech Stack:** Python 3.14・numpy・msgpack・lz4・pytest。MCP は FastMCP。上流の Key 番号の出典は MS-DIAL（`C:\Users\yuu18\MsdialWorkbench`、commit 45a531c）。

**Spec:**
- `docs/superpowers/specs/2026-11-06-unified-dataset-and-tools-design.md`（B。この計画は §4・§5・§7.1〜§7.3 を実装する。§6 定量は B-3、§7.4 設定と §8〜§10 公開面・削除・文書は B-2）
- `docs/superpowers/specs/2026-11-06-overhaul-overview-design.md`（全体。統合ブランチの方針）

## 設計との食い違い・設計が任せたことの決定（コードと上流で確かめたこと）

1. **同定の証拠の水準（§4.2 が「実装者が確かめて決める」とした欄）**: MS-DIAL の mzTab 書き出し（`src/MSDIAL5/MsdialCore/Export/MztabFormatExport.cs` の `ShouldWriteSmeLine`）を移植する。SME 相当（`feature_candidates`）になる条件は全部が真: (a) `IsMsmsAssigned`（`.arf` のそのスポットの試料のどれかで Key10 `MS2RawSpectrumID2CE` が空でない）、(b) 代表（`.arf2` Key56 の `representative`）の `Source` に Manual（64）が立っていない、(c) 代表が Text DB 由来でない（`Source & 16`、またはコンテナ Key2 `TextDbBasedMatchResults` に代表と同じ結果がある）、(d) 名前（Key12）が空でない、(e) `IsBlankFilteredByPostCurator`（Key58）が偽、(f) 名前に `no MS2` を含まない。上流にある「内部標準の正規化をしたのに内部標準が無い」の条件は `.mddata` のパラメータを読まないので移植しない（照合 Task 7 の既知の差に挙げる）。SML 相当（`feature_annotations`）は名前が `Unknown` / `null` / 空でない全スポット（SME 相当のスポットも含む。mzTab-M の SML と同じ）。
2. **特徴 ID**: `.arf` のピークの群には ID が無く位置で並ぶ（`arf/reader.py` の `_convert_arf_peak_group` は位置を `MasterAlignmentID` に入れている）。`.arf2` の Key0 `MasterAlignmentID` を `str()` したものを特徴 ID にし、`.arf` の群の数と `.arf2` のスポットの数が違えば読まない（Review Focus 1）。
3. **試料の順と名前**: 列は `.arf` の Key0 `FileID` の昇順。表示名は `.mddata` の今の名前（`FileID` で引く）を優先し、無ければ `.arf` の Key1 `FileName`。両者が違えば試料の情報に `arf_file_name` を残し、注記する（`docs/schema/AlignmentChromPeakFeature.md` の Key1 の注: アライン後の名前の変更で食い違いうる）。
4. **試料の種別**: `.mddata` の `AnalysisFileType`（0=Sample・1=Standard・2=QC・3=Blank。上流 `CommonEnums.cs`）を role に写す: `sample` / `standard` / `qc` / `blank`（`standard` は sample-manifest.v2 の role。生物試料として群の比較に入らない）。バッチは `AnalysisBatch`（Key11）、注入順は `AnalysisFileAnalyticalOrder`（Key4）、群は `AnalysisFileClass`（Key3）。`.mddata` が無ければ今の `detect_sample_roles`（名前から推す）に落とし、そのことを注記する。
5. **`.mddata` の Key 番号表が `docs/schema/` に無い**: Task 2 で上流 `src/MSDIAL5/MsdialCore/DataObj/AnalysisFile.cs`（`AnalysisFileBean`）から `docs/schema/AnalysisFileBean.md` を作ってから Key11 を読む（CLAUDE.md: 推測で直さない）。
6. **ハンドルの上限と対話の読み込み**: 読み込みのたびにハンドルが増える。上限（今の `LIPIDMIX_MAX_DATASET_HANDLES`、既定 8）に達したら今と同じ `DATASET_HANDLE_LIMIT` で止め、古いハンドルを黙って閉じない（閉じられた側の呼び出しが別のデータを読む事故を避ける。私の判断）。
7. **対話の呼び出しの並行**: 設計 §5 は「対話のセッションも同じ規則で」を実装者に任せた。`handle_dispatch` の包みで、`dataset_handle` を省いた呼び出しも直近のハンドルがあればそのロックの下でワーカースレッドに回す（同じデータセットへの対話と再解析の呼び出しが重ならない）。直近のハンドルが無ければ今どおりイベントループ上で動く。
8. **job / pipeline を解析フォルダのアダプタで読むときの出所の裏取り**: mzTab-M のアダプタの完了の門（`INCOMPLETE_ANALYSIS_JOB`・`exploratory_only`）と裏取り（`verified` / `legacy_unverified`）を、`loading.py` から関数に切り出して両方のアダプタで使う。解析フォルダのアダプタは、同じ job の mzTab-M の裏取りの結果をそのまま持つ（`.arf` は同じ run の出力として mzTab-M の隣から探すので、裏取りは run の単位）。
9. **pipeline を解析フォルダで読むとき**: 群は run が解いた試料対応表（`session_handoff._run_manifest_rows`）を使う。run の解析行列と差次的結果は載せない（D §6。v2 の行列は D で退役する）。`source="mztab"` は今と同じ（行列も差次的結果の復元も今どおり）で、再解析の既存の spec の数字を変えない。
10. **ライブラリの記憶の置き場所**: `DatasetState.library_record`（下の Interfaces）。パスはデータセットの中にだけ持ち、戻り値にはファイル名・種類・状態だけを出す（外部流出禁止）。
11. **`library="auto"` の意味の変更**: `dataset_handle_open` の `library="auto"` は「データセットが覚えたライブラリ」を開く（極性の設定から引く今の経路は、覚えたライブラリが無いときにも使わない。黙って違う版を使わないため）。極性の設定を下流から消すこと（`library_load` の `ion_mode`）と設定ファイルの一本化は B-2。

## Global Constraints

- コード・コメント・docstring・文書・コミットメッセージは日本語。
- テストは worktree の中から main ツリーの `.venv` で回す: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。
- コミットは pre-commit が全テストを回す（数分）。自動化した `git commit` はバックグラウンドで実行し、出力をファイルへ落として `FAILED` を引く。実行中にファイルを編集しない。`git stash` を使わない。作業は `overhaul/integration` の worktree（`.worktrees/overhaul-integration`）だけで行い、main を切り替えない。
- MessagePack の Key 番号は `docs/schema/*.md` から取る。表に無い Key を読むときは、先に上流の C# から表を足す。
- 研究室の参照データ（ライブラリ）のパス・ファイル名を追跡対象のテスト・文書に書かない。戻り値に出すのはファイル名まで（今の規則）。
- ツールの戻り値は `json_payload()`、全ツール `structured_output=False`、float は丸めて返す。
- fixture はテスト自身が tmp に作る（`tests/curation_fixtures.py` の `arf_row` `write_arf` `arf2_spot_raw` `write_arf2` `match_result` を使う）。実データのフォルダにテストを依存させない。
- 解析フォルダのアダプタが作る `DatasetState` は `detected_mask` を必ず持つ（`None` にしない）。gap-fill のセルは値を残し、検出状態だけを偽にする。
- 閉じた・存在しないハンドルを、対話のセッションや別のハンドルへ黙って落とさない（`DATASET_HANDLE_NOT_FOUND` か `missing_state`）。
- `source="mztab"` で読んだデータセットは、今の `dataset_load` / `dataset_handle_open` と同じ中身になる（analysis-spec.v1 などの再現のため。既存のテストがそのまま緑であること）。
- `metabolomix/core/{mcp_core,user_config,atomic_io,process_control}.py` に重い import を持ち込まない。`metabolomix/pipeline/{engine,worker,service,recovery,store}.py` から session を import しない。
- 新しい環境変数を足さない。

## Review Focus

1. **`.arf` と `.arf2` が別のバッチ・別のアラインメント**（スポット数が違う、または同じ数でも m/z が食い違う）→ 読まずに `ALIGNMENT_MISMATCH` で止める。位置で黙って重ねない。（Task 4 のテスト `test_mismatched_arf_and_arf2_refused`）
2. **アライン後に試料の名前を変えたフォルダ**（`.arf` の `FileName` と `.mddata` の名前が違う）→ `FileID` で結び、群とバッチを取り違えない。（Task 4 のテスト `test_renamed_sample_joined_by_file_id`）
3. **対話の直近のデータセットを `dataset_handle_close` で閉じた後**の `dataset_status` → `missing_state`。開いている別のハンドルを使わない。（Task 9 のテスト `test_closing_current_handle_does_not_fall_back`）
4. **覚えたライブラリのファイルが後から差し替わった**（sha256 が違う）→ 状態 `mismatch` で開かず、極性の設定にも落ちない。（Task 10 のテスト `test_library_sha_mismatch_not_opened`）
5. **名前が `Unknown` / `no MS2: …` / 空のスポット、Text DB 由来・手で直した代表** → SME 相当にしない（`w/o MS2: …` は上流の `ShouldWriteSmeLine` が除外しないので、MS/MS が付いていれば SME 相当。2026-11-06 の最終レビューで上流に合わせて直した。以前はここに `w/o MS2` を「SME 相当にしない」と書いていた）。`Unknown` と空は SML 相当にもしない。（Task 5 のテスト `test_sme_rule_ported_from_mztab_export`）

---

### Task 1: `.arf2` のスポットの欄（同定の写しに要る欄）

**Files:**
- Create: `metabolomix/arf2/spot_fields.py`
- Test: `tests/test_arf2_spot_fields.py`

**Interfaces:**
- Consumes: `metabolomix.arf2.reader.load_raw_spots(path) -> list[list]`、`_extract_times`、`metabolomix.arf2.match_results.representative(container)` / `_candidates(container)` / `decode_match_result(raw)`
- Produces:
  - `spot_fields(raw_spot: list) -> dict`。キー: `master_id: int`（Key0）、`name: str | None`（Key12）、`formula: str | None`（Key13 の [0]）、`ontology`（Key14）、`smiles`（Key15）、`inchikey`（Key16）、`adduct: str | None`（Key54 の [2]）、`mz: float | None`（Key5）、`rt: float | None`（Key4 の ChromXs の RT、分）、`ion_mode: str`（Key11 を `Positive` / `Negative` / `Both` / `Unknown`）、`comment: str | None`（Key25）、`internal_standard_id: int | None`（Key7）、`blank_filtered: bool`（Key58）、`fill_percentage: float | None`（Key49）、`representative: dict | None`、`is_text_db_representative: bool`、`is_manually_modified: bool`、`n_candidates: int`。空文字列は `None` にそろえる。
  - `load_spot_fields(arf2_path) -> list[dict]`（ファイル順）

- [ ] **Step 1: 失敗するテストを書く**

```python
from tests.curation_fixtures import arf2_spot_raw, match_result, write_arf2
from metabolomix.arf2.spot_fields import load_spot_fields, spot_fields

def test_spot_fields_reads_identity_keys(tmp_path):
    raw = arf2_spot_raw(spot_id=7, name="PC 34:1", mz=760.5851, rt=12.0,
                        adduct="[M+H]+", matches=[match_result()])
    raw[25] = "checked"; raw[7] = 3; raw[58] = False
    f = spot_fields(raw)
    assert f["master_id"] == 7 and f["name"] == "PC 34:1"
    assert f["adduct"] == "[M+H]+" and f["ontology"] == "PC"
    assert f["comment"] == "checked" and f["internal_standard_id"] == 3
    assert f["rt"] == 12.0 and f["ion_mode"] == "Positive"
    assert f["representative"]["name"] == "PC 34:1"
    assert f["is_text_db_representative"] is False and f["is_manually_modified"] is False

def test_text_db_and_manual_representative_flags():
    text = arf2_spot_raw(spot_id=1, matches=[match_result({26: 16})])
    assert spot_fields(text)["is_text_db_representative"] is True
    in_list = arf2_spot_raw(spot_id=2, matches=[match_result()])
    in_list[56][2] = [match_result()]          # コンテナ Key2 TextDbBasedMatchResults
    assert spot_fields(in_list)["is_text_db_representative"] is True
    manual = arf2_spot_raw(spot_id=3, matches=[match_result({26: 4 | 64})])
    assert spot_fields(manual)["is_manually_modified"] is True

def test_load_spot_fields_keeps_file_order(tmp_path):
    path = write_arf2(tmp_path / "A.arf2", [arf2_spot_raw(spot_id=i) for i in (0, 1, 2)])
    assert [f["master_id"] for f in load_spot_fields(path)] == [0, 1, 2]
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_arf2_spot_fields.py -q`
Expected: FAIL（`ModuleNotFoundError: metabolomix.arf2.spot_fields`）

- [ ] **Step 3: `spot_fields` と `load_spot_fields` を実装する**

Key の出典は `docs/schema/AlignmentSpotProperty.md` と `docs/schema/MsScanMatchResult.md`（コンテナは Key0 結果の一覧・Key1 MSP の辞書・Key2 Text DB の一覧）。`is_text_db_representative` は上流 `MsScanMatchResultContainer.IsTextDbBasedRepresentative` の移植（`Source & 16`、または Key2 の一覧に代表と同じ復号結果がある）。`is_manually_modified` は代表の `is_manually_modified`（`match_results.decode_match_result` が既に出す）。

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_arf2_spot_fields.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add metabolomix/arf2/spot_fields.py tests/test_arf2_spot_fields.py
git commit -m "feat(arf2): 同定の写しに要るスポットの欄を読む spot_fields を足す"
```

---

### Task 2: `.mddata` の試料の欄（バッチ）と Key 番号表

**Files:**
- Create: `docs/schema/AnalysisFileBean.md`
- Modify: `metabolomix/msdial/classes.py`（`parse_analysis_file_classes` と Key 定数）
- Create: `tests/msdial_folder_fixtures.py`（この Task では `.mddata` を書く部分だけ）
- Modify: `tests/test_msdial_classes.py`（`pack_msdial` と `analysis_file_row` を fixture から import する）

**Interfaces:**
- Produces:
  - `parse_analysis_file_classes(mddata_path) -> list[dict]` の各行に `analysis_batch: int | None`（Key11）を足す。既存のキーは変えない。
  - `tests/msdial_folder_fixtures.py`: `pack_msdial(value) -> bytes`、`analysis_file_row(file_id: int, name: str, class_id: str, *, file_type: int = 0, order: int | None = None, batch: int = 1) -> list`、`write_mddata(path, rows: list[list]) -> Path`

- [ ] **Step 1: Key 番号表を作る**

上流 `src/MSDIAL5/MsdialCore/DataObj/AnalysisFile.cs`（commit 45a531c、2026-09-02）の `[Key(0)]`〜`[Key(18)]` を、`docs/schema/AlignmentSpotProperty.md` と同じ書式で表にする。`AnalysisFileType` の列挙（`src/Common/CommonStandard/Enum/CommonEnums.cs`: 0=Sample・1=Standard・2=QC・3=Blank）も書く。

- [ ] **Step 2: 失敗するテストを書く**（`tests/test_msdial_classes.py` に足す）

```python
def test_parse_reads_analysis_batch(tmp_path):
    path = write_mddata(tmp_path / "d.mddata",
                        [analysis_file_row(0, "s1", "A", batch=2, file_type=2, order=5)])
    row = parse_analysis_file_classes(path)[0]
    assert row["analysis_batch"] == 2
    assert row["analysis_file_type"] == 2 and row["analytical_order"] == 5
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_msdial_classes.py -q`
Expected: FAIL（`ImportError` か `KeyError: 'analysis_batch'`）

- [ ] **Step 4: 実装する**

`classes.py` に `ANALYSIS_BATCH_INDEX = 11` を足し、`parse_analysis_file_classes` の行に `analysis_batch` を足す。fixture の `analysis_file_row` は `test_msdial_classes.py` の今の関数を移し、Key2・Key4・Key11 を引数で置けるようにする（既定は今の値）。

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_msdial_classes.py -q`
Expected: PASS

- [ ] **Step 6: コミット**

```bash
git add docs/schema/AnalysisFileBean.md metabolomix/msdial/classes.py tests/msdial_folder_fixtures.py tests/test_msdial_classes.py
git commit -m "feat(msdial): .mddata の AnalysisBatch を読み、AnalysisFileBean の Key 番号表を置く"
```

---

### Task 3: 解析フォルダのファイルの選択

**Files:**
- Create: `metabolomix/dataset/__init__.py`（空）、`metabolomix/dataset/alignment_files.py`
- Modify: `CLAUDE.md`（構成の木に `metabolomix/dataset/` の 1 行）
- Test: `tests/test_alignment_files.py`

**Interfaces:**
- Consumes: `metabolomix.core.path_resolvers._select_latest_batch` / `_pick_latest` / `_batch_key`、`metabolomix.msdial.alignment_spot.alignment_sibling(path, suffix)`、`metabolomix.msdial.classes.resolve_mddata_path(path)`、`metabolomix.msdial.tags.alignment_tag_path(arf2_path)`
- Produces:
  - `@dataclass(frozen=True) class AlignmentSelection`: `arf: Path`、`arf2: Path`、`dcl: Path | None`、`eic: Path | None`、`tags: Path | None`、`mddata: Path | None`、`loaded_libraries: tuple[Path, ...]`（同じフォルダの `*_Loaded.*.dbs`）、`not_chosen: tuple[str, ...]`（選ばなかった `.arf` のファイル名）
  - `select_alignment(path: str | Path) -> AlignmentSelection`。`path` はフォルダか `.arf` / `.arf2` のファイル。
  - エラー: `DomainError("ALIGNMENT_NOT_FOUND", ...)`（`.arf` が無い）、`DomainError("ALIGNMENT_INCOMPLETE", ...)`（同じバッチの `.arf2` が無い）

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_prefers_peak_properties_and_latest_batch(tmp_path):
    for stem in ("AlignmentResult_2026_01_01_00_00_00", "AlignmentResult_2026_02_01_00_00_00"):
        for suffix in ("_PeakProperties.arf", "_DriftSpots.arf", ".arf2", ".dcl"):
            (tmp_path / f"{stem}{suffix}").write_bytes(b"x")
    sel = select_alignment(tmp_path)
    assert sel.arf.name == "AlignmentResult_2026_02_01_00_00_00_PeakProperties.arf"
    assert sel.arf2.name == "AlignmentResult_2026_02_01_00_00_00.arf2"
    assert sel.dcl.name == "AlignmentResult_2026_02_01_00_00_00.dcl"
    assert "AlignmentResult_2026_01_01_00_00_00_PeakProperties.arf" in sel.not_chosen

def test_arf_without_arf2_is_incomplete(tmp_path):
    (tmp_path / "AlignmentResult_2026_02_01_00_00_00_PeakProperties.arf").write_bytes(b"x")
    with pytest.raises(DomainError) as exc:
        select_alignment(tmp_path)
    assert exc.value.code == "ALIGNMENT_INCOMPLETE"

def test_no_arf_is_not_found(tmp_path):
    with pytest.raises(DomainError) as exc:
        select_alignment(tmp_path)
    assert exc.value.code == "ALIGNMENT_NOT_FOUND"

def test_console_layout_alignresult_names(tmp_path):
    for suffix in (".arf", ".arf2", ".dcl"):
        (tmp_path / f"AlignResult-2026_03_01_10_00_00{suffix}").write_bytes(b"x")
    sel = select_alignment(tmp_path / "AlignResult-2026_03_01_10_00_00.arf2")
    assert sel.arf.name == "AlignResult-2026_03_01_10_00_00.arf"
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_alignment_files.py -q`
Expected: FAIL（`ModuleNotFoundError`）

- [ ] **Step 3: `select_alignment` を実装する**

フォルダなら直下の `*.arf` を `_select_latest_batch` で最新バッチへ絞り、`PeakProperties` を含む名前を優先して `_pick_latest`。`.arf2` を渡されたら同じ幹の `.arf`（`_PeakProperties.arf` → `.arf` の順）を探す。`.arf2` `.dcl` `.EIC.aef` は `alignment_sibling`。サブフォルダは探さない（別の解析を掴まない）。

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_alignment_files.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add metabolomix/dataset tests/test_alignment_files.py CLAUDE.md
git commit -m "feat(dataset): 解析フォルダから同じバッチのアラインメントの一式を選ぶ"
```

---

### Task 4: 解析フォルダのアダプタ（値・検出状態・試料・極性）

**Files:**
- Create: `metabolomix/dataset/msdial_adapter.py`
- Modify: `metabolomix/mztab/dataset_state.py`（欄の追加）、`metabolomix/analysis/sample_manifest.py`（`_default_row`）、`metabolomix/mztab/polarity.py`（`dataset_polarity`）
- Modify: `tests/msdial_folder_fixtures.py`（フォルダ一式を書く）
- Test: `tests/test_msdial_adapter.py`

**Interfaces:**
- Consumes: Task 1 `load_spot_fields`、Task 2 `parse_analysis_file_classes`、Task 3 `AlignmentSelection`、`metabolomix.arf.reader.deserialize` / `alignment_feature_row`、`metabolomix.mztab.evidence.compare_feature_mz`、`sample_manifest.resolve_metadata` / `apply_metadata`
- Produces:
  - `load_msdial_folder(selection: AlignmentSelection, *, measure: str = "peak_height") -> DatasetState`。`measure` は `peak_height`（`.arf` Key18）か `peak_area_above_zero`（Key20）。それ以外は `DomainError("DATASET_BAD_REQUEST")`。
  - `DatasetState` に足す欄（`__init__` で既定値を置く）: `msdial_samples: list[dict] | None`（試料ごとに `file_id` `name` `arf_file_name` `group` `role` `file_type` `batch` `injection_order`）、`ion_mode_majority: str | None`（`positive` / `negative`）、`artifact_sha256: dict[str, str]`、`library_record: dict | None`（Task 6）、`feature_tags: dict[str, list[str]]`・`curation_flags: dict | None`（Task 5）、`load_notes: list[str]`。
  - 埋める既存の欄: `source_format="msdial_folder"`、`quantification_measure`、`feature_ids`、`sample_names`、`sample_assay_ids`（`file_id` の文字列）、`feature_matrix`（特徴×試料）、`detected_mask`、`per_sample_mz`、`feature_metadata`（この Task では `name` `mz` `rt` `ontology` `adduct` `formula` `comment` `fill_percentage` `internal_standard_id` `msms_assigned`）、`source_files`（`.arf` `.arf2` `.mddata` の sha256）、`artifact_paths`（`{"dcl": [...], "eic": [...], "tags": [...]}`）、`roles`、`sample_meta`、`sample_metadata_rows`（`apply_metadata` 経由）。
  - `tests/msdial_folder_fixtures.py`: `write_msdial_folder(folder: Path, *, spots: list[dict], samples: list[dict], with_mddata: bool = True, stem: str = "AlignmentResult_2026_02_01_00_00_00") -> dict`（`spots[i]` は `{"name","mz","rt","ontology","adduct","heights":[...],"gap_filled":[...],"msms":[...],"matches":[...]}`、`samples[j]` は `{"file_id","name","class_id","file_type","order","batch"}`。戻り値は書いたパスの辞書）

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_values_detection_and_ids(tmp_path):
    paths = write_msdial_folder(tmp_path, spots=[
        {"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0],
         "gap_filled": [False, True, False]},
        {"name": "Unknown", "mz": 500.1, "heights": [10.0, 20.0, 30.0],
         "gap_filled": [False, False, False]}],
        samples=[{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(3)])
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.source_format == "msdial_folder" and ds.feature_ids == ["0", "1"]
    assert ds.sample_names == ["s0", "s1", "s2"]
    assert ds.feature_matrix[0].tolist() == [100.0, 0.0, 50.0]
    assert ds.detected_mask[0].tolist() == [True, False, True]
    assert "s1" not in ds.per_sample_mz["0"]
    assert ds.quantification_measure == "peak_height"

def test_sample_roles_batches_and_groups_from_mddata(tmp_path):
    write_msdial_folder(tmp_path, spots=[...one spot...], samples=[
        {"file_id": 0, "name": "QC1", "class_id": "QC", "file_type": 2, "order": 1, "batch": 1},
        {"file_id": 1, "name": "x", "class_id": "A", "file_type": 0, "order": 2, "batch": 2},
        {"file_id": 2, "name": "std", "class_id": "S", "file_type": 1, "order": 3, "batch": 2},
        {"file_id": 3, "name": "blk", "class_id": "B", "file_type": 3, "order": 4, "batch": 2}])
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.roles == {"QC1": "qc", "x": "sample", "std": "standard", "blk": "blank"}
    rows = {r["sample_id"]: r for r in ds.sample_metadata_rows}
    assert rows["x"]["group"] == "A" and rows["x"]["batch"] == "2" and rows["x"]["injection_order"] == 2
    assert rows["x"]["provenance"]["group"]["source"] == "msdial_mddata"

def test_without_mddata_falls_back_to_name_roles(tmp_path):
    write_msdial_folder(tmp_path, spots=[...], samples=[...], with_mddata=False)
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert any("mddata" in note for note in ds.load_notes)

def test_mismatched_arf_and_arf2_refused(tmp_path):          # Review Focus 1
    paths = write_msdial_folder(tmp_path, spots=[...2 spots...], samples=[...])
    write_arf2(paths["arf2"], [arf2_spot_raw(spot_id=0)])   # スポットを 1 件に減らす
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ALIGNMENT_MISMATCH"

def test_renamed_sample_joined_by_file_id(tmp_path):         # Review Focus 2
    # .arf の FileName は "old_b"、.mddata の今の名前は "s1"（FileID 1）
    ...
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.sample_names[1] == "s1"
    assert ds.msdial_samples[1]["arf_file_name"] == "old_b"
    assert {r["sample_id"]: r["group"] for r in ds.sample_metadata_rows}["s1"] == "B"

def test_area_measure_and_bad_measure(tmp_path): ...       # peak_area_above_zero は Key20、"foo" は DATASET_BAD_REQUEST

def test_polarity_from_arf2_ion_mode(tmp_path):
    ...  # 全スポット ion_mode=1 → dataset_polarity(ds) == ("negative", "arf2_ion_mode", [])
```

（`...` の箇所は `write_msdial_folder` の引数を同じ形で書く。各テストの中身はこの見出しの主張だけを確かめる。）

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_msdial_adapter.py -q`
Expected: FAIL

- [ ] **Step 3: アダプタを実装する**

- `.arf` は `arf_reader.deserialize` で群を読み、各セルを `alignment_feature_row` で取る（Key を直に書かない）。群の数 ≠ `.arf2` のスポット数なら `ALIGNMENT_MISMATCH`。数が同じでも `compare_feature_mz`（群の代表 `MassCenter` と `.arf2` の Key5）が不合格なら `ALIGNMENT_MISMATCH`（理由に最悪の差を入れる）。
- 試料の列は全群で現れた `file_id` の昇順。あるセルが欠けた群は値 `nan`・検出 False。
- `.mddata` の行は `file_id` で引く。`.mddata` に有って `.arf` に無い試料は使わない。
- 試料の行は `sample_manifest._default_row` を直して作る: `ds.msdial_samples` があればその `role` `group` `batch`（文字列）`injection_order` を使い、provenance は `{"source": "msdial_mddata", "confidence": "declared"}`。無ければ今の推定のまま。最後に `apply_metadata(ds, resolve_metadata(ds, None))`。
- `dataset_polarity`: `job_path` の宣言を最優先にする今の順は変えず、宣言が無く `ds.source_format == "msdial_folder"` なら `ds.ion_mode_majority` を出所 `arf2_ion_mode` で返す（mzTab の付加体の多数決は見ない）。`ion_mode_majority` は `.arf2` の `ion_mode` の多数（`Both` / `Unknown` は数えない、陽陰が同数なら None）。
- `source_files` と `artifact_sha256` の sha256 は `dataset_state._sha256` を使う。

- [ ] **Step 4: 通ることを確かめる**（既存の sample_manifest と polarity のテストも）

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_msdial_adapter.py tests/test_sample_manifest.py tests/test_dataset_handle_tools.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add metabolomix/dataset/msdial_adapter.py metabolomix/mztab/dataset_state.py metabolomix/analysis/sample_manifest.py metabolomix/mztab/polarity.py tests/msdial_folder_fixtures.py tests/test_msdial_adapter.py
git commit -m "feat(dataset): 解析フォルダから値・検出状態・試料・極性を持つ DatasetState を作る"
```

---

### Task 5: アダプタの同定・標準化・タグ・キュレーション

**Files:**
- Create: `metabolomix/dataset/msdial_identity.py`
- Modify: `metabolomix/dataset/msdial_adapter.py`（呼び出し）
- Test: `tests/test_msdial_identity.py`

**Interfaces:**
- Consumes: Task 1 の `spot_fields` の辞書、Task 4 の `msms_assigned`、`metabolomix.mztab.dataset_state._all_candidates` と `_build_feature_annotations`（作る辞書のキーの正準。読むだけ）、`metabolomix.core.tool_helpers._identity_tables()`、`metabolomix.msdial.lipid_identity.build_identity_block`、`metabolomix.msdial.tags.parse_tag_file` / `alignment_tag_path`、`metabolomix.curation.apply.flags_for_arf2`
- Produces:
  - `is_sme_like(fields: dict, *, msms_assigned: bool) -> bool`（「食い違い 1」の (a)〜(f)）
  - `sme_candidate(fields: dict, fid: str) -> dict`（キーの集合は mzTab の `_all_candidates` の 1 件と同じ。`rank=1`、`sme_id=f"arf2:{fid}"`、`confidence_value` は代表の `total_score`）
  - `sml_annotation(fields: dict) -> dict | None`（キーの集合は `_build_feature_annotations` の 1 件と同じ。曖昧さは持たないので `ambiguous=False`）
  - `identity_columns(fields: dict, *, msms_assigned: bool, tables) -> dict` → `{"normalized", "refmet", "lipid_maps_category", "msi_level"}`（`build_identity_block` を `mass_error_band="UNKNOWN"`、`adduct_band="UNKNOWN"` で。同じ `(name, ontology, msms_assigned)` は 1 回だけ計算する）
  - アダプタが埋める: `feature_candidates`、`feature_annotations`、`feature_metadata[fid]` の `inchikey` `inchikey_source="arf2"` `smiles` `inchi=None` `annotator_id` `annotation_source`（代表の Source を `msp` / `text_db` / `generated_lipid` / `manual` / `none`）`spectrum_matched`（代表の `is_spectrum_match`。代表が無ければ None）`total_score`（代表の `total_score`、丸めない）`identity`（`identity_columns` の辞書）、`feature_tags`（`{fid: [ラベル, ...]}`。タグファイルが無ければ空）、`curation_flags`（`flags_for_arf2` の戻り値。`FlagFileError` は読み込みを止めず `load_notes` に 1 行）

- [ ] **Step 1: 失敗するテストを書く**

```python
@pytest.mark.parametrize("name,overrides,msms,expected", [
    ("PC 34:1", {}, True, True),
    ("PC 34:1", {}, False, False),                 # (a) MS/MS 未取得
    ("PC 34:1", {26: 4 | 64}, True, False),        # (b) 手で直した代表
    ("PC 34:1", {26: 16}, True, False),            # (c) Text DB
    ("", {}, True, False),                         # (d) 名前が空
    ("no MS2: PC 34:1", {}, True, False),          # (f)
])
def test_sme_rule_ported_from_mztab_export(name, overrides, msms, expected):   # Review Focus 5
    raw = arf2_spot_raw(spot_id=0, name=name, matches=[match_result(overrides)])
    assert is_sme_like(spot_fields(raw), msms_assigned=msms) is expected

def test_blank_filtered_spot_is_not_sme():
    raw = arf2_spot_raw(spot_id=0, matches=[match_result()]); raw[58] = True
    assert is_sme_like(spot_fields(raw), msms_assigned=True) is False

def test_candidate_and_annotation_keys_match_mztab_path(tmp_path):
    ds_mz = <reanalysis_fixtures.write_mztab で作った mzTab を load_dataset_state で読む>
    fid = next(iter(ds_mz.feature_candidates))
    raw = arf2_spot_raw(spot_id=0, matches=[match_result()])
    assert set(sme_candidate(spot_fields(raw), "0")) == set(ds_mz.feature_candidates[fid][0])
    fid_a = next(iter(ds_mz.feature_annotations))
    assert set(sml_annotation(spot_fields(raw))) == set(ds_mz.feature_annotations[fid_a])

def test_unknown_has_no_annotation():
    assert sml_annotation(spot_fields(arf2_spot_raw(spot_id=0, name="Unknown"))) is None

def test_adapter_fills_identity_tags_and_flags(tmp_path):
    write_msdial_folder(tmp_path, ...)       # spot 0 は MS/MS あり・参照一致、spot 1 は Unknown
    <_tags.xml に spot 0 の Confirmed を書く（parse_tag_file が読む形）>
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_candidates["0"][0]["chemical_name"] == "PC 34:1"
    assert "1" not in ds.feature_candidates and "1" not in ds.feature_annotations
    assert ds.feature_metadata["0"]["identity"]["msi_level"] in (2, 3)
    assert ds.feature_tags["0"] == ["Confirmed"]
    assert ds.curation_flags["n"] == 0
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_msdial_identity.py -q`
Expected: FAIL

- [ ] **Step 3: 実装する**

`sme_candidate` の値は代表の照合結果（`name` → `chemical_name`、`inchikey`、`total_score`）とスポット（`adduct`、`formula`、`smiles`、`mz` → `exp_mass_to_charge`、`ontology`、`spectrum_matched` は代表の `is_spectrum_match`）から取り、mzTab にしかない欄（`database_identifier` など）は `None`。

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_msdial_identity.py tests/test_msdial_adapter.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add metabolomix/dataset/msdial_identity.py metabolomix/dataset/msdial_adapter.py tests/test_msdial_identity.py
git commit -m "feat(dataset): 解析フォルダの同定を mzTab-M と同じ規則で SME/SML 相当に写し、標準化・タグ・キュレーションを載せる"
```

---

### Task 6: ライブラリの記憶

**Files:**
- Create: `metabolomix/dataset/library_memory.py`
- Modify: `metabolomix/mztab/loading.py`（実効メソッドから MSP の場所を引く部分を関数に切り出す）
- Modify: `docs/output_format/library.md`（記憶の状態と LBM の印の意味）
- Test: `tests/test_library_memory.py`

**Interfaces:**
- Consumes: Task 3 `AlignmentSelection.loaded_libraries`、`loading._msp_from_annotator_table` / `read_method_keys`、`pipeline.store.load_run(root)` の `record["inputs"]` の `msp`（`{"path","file","sha256","source"}`）・`lbm`・`library_mode`
- Produces:
  - `loading.effective_method_msp(job) -> Path | None`（`_msp_table_for_job` の MSP を探す部分。`_msp_table_for_job` はこれを呼ぶように直す）
  - `LibraryRecord`（dict）: `status`（`known` / `unknown` / `mismatch` / `lbm_only`）、`kind`（`dbs` / `msp` / `lbm` / None）、`file`（ファイル名）、`path`（データセットの中だけで持つ。戻り値に出さない）、`sha256`、`origin`（`loaded_dbs` / `effective_method` / `run_inputs` / None）、`candidates`（ファイル名の一覧）、`lbm_also: bool`（`msp_then_lbm`）
  - `remember_from_folder(selection) -> dict`、`remember_from_job(job_path: Path) -> dict`、`remember_from_run(pipeline_root: Path) -> dict`、`unknown_record() -> dict`（mzTab-M を直接読んだとき。`status="unknown"`、ほかは None / 空）
  - `verify_record(record: dict) -> dict`（今のファイルの sha256 を測り直し、違えば `status="mismatch"` の写しを返す。ファイルが無ければ `mismatch`）
  - `public_view(record: dict | None) -> dict | None`（`path` を外した写し）
  - `is_lbm_annotator(annotator_id: str | None) -> bool | None`（Step 3 で決める。決められなければ常に None を返す関数として置く）

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_folder_with_one_loaded_dbs_is_known(tmp_path):
    write_msdial_folder(tmp_path, ...)
    write_dbs(tmp_path / "Proj_Loaded.msp2.dbs", [...])         # tests/dbs_fixture.py
    rec = remember_from_folder(select_alignment(tmp_path))
    assert rec["status"] == "known" and rec["kind"] == "dbs" and rec["origin"] == "loaded_dbs"
    assert "path" not in public_view(rec)

def test_folder_with_two_loaded_dbs_is_unknown_with_candidates(tmp_path): ...
def test_folder_with_only_lbm_dbs_is_lbm_only(tmp_path): ...    # *_Loaded.lbm2.dbs だけ
def test_job_records_effective_method_msp(tmp_path): ...        # 実効メソッドが MSP を宣言 → kind "msp", origin "effective_method"
def test_run_pinned_sha_mismatch(tmp_path):
    ...  # record["inputs"]["msp"]["sha256"] と今のファイルが違う → status "mismatch"
def test_msp_then_lbm_records_msp_with_lbm_also(tmp_path): ...
def test_mztab_direct_is_unknown():
    assert unknown_record()["status"] == "unknown"               # Task 8 が mzTab-M の直接読みで使う
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_library_memory.py -q`
Expected: FAIL

- [ ] **Step 3: LBM の印の規則を決める（調べる手順）**

`library_mode="msp_then_lbm"` の reanalysis の run（`C:\Users\yuu18\datasets\reanalysis\` の下、`pipeline-run.json` の `inputs.library_mode` で探す）の `.arf2` を `load_spot_fields` で読み、代表の `annotator_id` と `source` を MSP 由来と LBM 由来で比べる。上流（`MsdialWorkbench` の `AnnotatorID` を付けるところ）で接頭辞などの規則が確かめられれば、その規則で `is_lbm_annotator` を実装し、テストに規則の値（ライブラリのファイル名を含まない部分）を書く。確かめられなければ `is_lbm_annotator` は常に None を返し、`docs/output_format/library.md` に「LBM が代表になった特徴の印は付けられない（理由）」を書く。どちらにしたかを `docs/HISTRY.md`（main ツリー）に書く。

- [ ] **Step 4: 実装し、通ることを確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_library_memory.py tests/test_dataset_loading_service.py tests/test_dataset_state.py -q`
Expected: PASS（`_msp_table_for_job` の既存のテストも緑）

- [ ] **Step 5: コミット**

```bash
git add metabolomix/dataset/library_memory.py metabolomix/mztab/loading.py docs/output_format/library.md tests/test_library_memory.py
git commit -m "feat(dataset): データセットが同定に使ったライブラリを覚える"
```

---

### Task 7: 2 つのアダプタの照合

**Files:**
- Create: `metabolomix/dataset/agreement.py`
- Modify: `tests/msdial_folder_fixtures.py`（同じアラインメントの mzTab-M を書く）
- Test: `tests/test_adapter_agreement.py`

**Interfaces:**
- Consumes: Task 4・5 の `DatasetState`（`source_format="msdial_folder"`）、mzTab-M のアダプタの `DatasetState`、`evidence.MZ_TOLERANCE`
- Produces:
  - `compare_states(folder_ds, mztab_ds, *, rel_tol: float = 1e-4) -> dict`。戻り値: `{"ok": bool, "checks": {"feature_ids": {"only_folder": [...≤5], "only_mztab": [...≤5], "n_common": int}, "samples": {"same_set": bool, "same_order": bool}, "values": {"n_compared": int, "n_over": int, "worst_rel": float}, "detected": {"n_compared": int, "n_diff": int}, "names": {"n_diff": int, "examples": [...≤5]}, "evidence": {"n_diff": int, "examples": [...≤5]}}, "known_differences": [str]}`。`ok` は特徴 ID と試料の集合が一致し、`values.n_over`・`detected.n_diff`・`names.n_diff`・`evidence.n_diff` がすべて 0。
  - `known_differences` には「食い違い 1」の移植しない条件（内部標準の正規化）を固定文で入れる。
  - `tests/msdial_folder_fixtures.py`: `write_matching_mztab(path: Path, folder_spec: dict) -> Path`（`write_msdial_folder` に渡したのと同じ `spots` / `samples` から、SMF_ID・assay 名・値・SME/SML を mzTab-M に書く。`tests/reanalysis_fixtures.write_mztab` の書式に合わせる）

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_same_alignment_agrees(tmp_path):
    spec = {"spots": [...3 spots: SME あり / SML だけ / Unknown...], "samples": [...4 samples...]}
    write_msdial_folder(tmp_path, **spec)
    mz = write_matching_mztab(tmp_path / "out.mzTab", spec)
    report = compare_states(load_msdial_folder(select_alignment(tmp_path)),
                            load_dataset_state(mztab_path=mz))
    assert report["ok"] is True, report

def test_value_difference_is_reported(tmp_path):
    ...  # mzTab 側の 1 セルだけ 10% 違う値で書く
    assert report["ok"] is False and report["checks"]["values"]["n_over"] == 1

def test_evidence_level_difference_is_reported(tmp_path):
    ...  # mzTab 側で SME を書かない → evidence.n_diff == 1
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_adapter_agreement.py -q`
Expected: FAIL

- [ ] **Step 3: `compare_states` を実装する**

値は共通の特徴 × 共通の試料で `|a-b| / max(|a|,|b|,1)` を `rel_tol` と比べる（両方 `nan` は一致）。検出状態は mzTab 側の `detected_mask` が `None` なら比べず `n_compared=0`。名前は SME の最上位の `chemical_name`、無ければ SML の `name`。証拠の水準は `sme` / `sml` / `none`（`feature_search._window_identity` と同じ決め方）。

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_adapter_agreement.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add metabolomix/dataset/agreement.py tests/msdial_folder_fixtures.py tests/test_adapter_agreement.py
git commit -m "feat(dataset): 解析フォルダと mzTab-M の 2 つのアダプタの DatasetState を突き合わせる"
```

---

### Task 8: 入口 `open_dataset` と既存ツールへの配線

**Files:**
- Create: `metabolomix/dataset/opening.py`
- Modify: `metabolomix/mztab/loading.py`（job の門の切り出し）、`metabolomix/tools/handle_tools.py`（`_signature` `_load` を `opening` へ移し、引数を足す）、`metabolomix/tools/mztab_tools.py`（`dataset_load` の引数）、`USAGE.md`（2 ツールの引数）
- Test: `tests/test_dataset_opening.py`

**Interfaces:**
- Consumes: Task 3〜6、`loading.load_dataset_state`、`session_handoff.load_run_dataset` / `_run_manifest_rows`、`mztab.evidence.arf_candidates(mztab_path)`、`handle_tools._pipeline_arf2`、`handles.registry.cached_base` / `put_base`
- Produces:
  - `loading.job_gate(job_path: Path, *, allow_incomplete: bool) -> JobGate`。`JobGate`（dataclass）: `job`、`mztab_path: Path`、`exploratory_only: bool`、`verified: bool`、`reason: str | None`。`_load_from_job` はこれを使うように直す（中身の規則は変えない）。
  - `@dataclass class OpenedDataset`: `dataset: DatasetState`、`kind: str`（`directory` / `mztab_path` / `job_path` / `pipeline_path`）、`path: str`、`adapter: str`（`msdial_folder` / `mztab`）、`selection: AlignmentSelection | None`、`notes: list[str]`
  - `open_dataset(*, directory=None, mztab_path=None, job_path=None, pipeline_path=None, source: str = "auto", measure: str = "peak_height", allow_incomplete: bool = False) -> OpenedDataset`
    - 4 つの入口のちょうど 1 つ（0 個・2 個以上は `DATASET_BAD_REQUEST`）。`source` は `auto` / `folder` / `mztab`。
    - `directory` → 解析フォルダ（`source="mztab"` は `DATASET_BAD_REQUEST`）。`mztab_path` → mzTab-M（`source="folder"` は `DATASET_BAD_REQUEST`）。`job_path` / `pipeline_path` → `auto` と `folder` は解析フォルダ、`mztab` は今の経路そのまま。
    - job / pipeline を解析フォルダで読むとき: `job_gate` で mzTab-M を決め、その隣の `.arf` を `arf_candidates` の先頭で探し（pipeline は無ければ `_pipeline_arf2`）、`select_alignment`。見つからなければ `ALIGNMENT_NOT_FOUND`（mzTab-M へ黙って落とさない。理由に `source="mztab"` を案内する）。`exploratory_only`・`source_verification`・`job_path` は門の結果を写す。pipeline は `_run_manifest_rows` があれば `apply_metadata(ds, resolve_metadata(ds, rows))`。
    - `library_record` は入口ごとに Task 6 の関数で埋める（mzTab-M を直接読んだときは `unknown_record()`。`source="mztab"` で job / pipeline を読んだときは `remember_from_job` / `remember_from_run`）。
    - `dataset_handle` を受けるツールに引数を足すこと（設計 §5 の「全ツールが省略可能な `dataset_handle` を持つ」）は B-2 で行う。この Task は入口の 2 ツールだけを変える。
    - 基のキャッシュ: キーは `f"{kind}:{adapter}:{measure}:{絶対パス}"`、指紋は今の `_signature`。返す `dataset` は常に `copy.deepcopy`。
  - `dataset_load(mztab_path=None, job_path=None, pipeline_path=None, allow_incomplete=False, directory=None, source="auto", measure="peak_height")` と `dataset_handle_open(..., directory=None, source="auto", measure="peak_height")`（既存の引数と並びは変えず、末尾に足す）。戻り値に `adapter`・選んだ `.arf` のファイル名・選ばなかった候補・`library_record` の `public_view` を足す。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_directory_opens_folder_adapter(tmp_path): ...       # adapter == "msdial_folder"
def test_exactly_one_entry_required(): ...                    # 0 個・2 個 → DATASET_BAD_REQUEST
def test_source_conflicts_rejected(tmp_path): ...             # directory+mztab / mztab_path+folder
def test_job_defaults_to_folder_and_mztab_on_request(tmp_path):
    ...  # job の run の隣に .arf 一式と mzTab-M を置く（tests/metabolomics_fixtures か console のテストの job の作り方に合わせる）
    assert open_dataset(job_path=job).adapter == "msdial_folder"
    assert open_dataset(job_path=job, source="mztab").adapter == "mztab"
def test_incomplete_job_gate_shared(tmp_path):
    ...  # status="partial" の job: 既定は INCOMPLETE_ANALYSIS_JOB、allow_incomplete=True なら両アダプタとも exploratory_only
def test_job_without_arf_does_not_fall_back_to_mztab(tmp_path): ...   # ALIGNMENT_NOT_FOUND
def test_mztab_source_matches_current_loader(tmp_path):
    ...  # source="mztab" の feature_matrix・feature_ids・sample_names が load_dataset_state と同じ
def test_base_cache_returns_independent_copies(tmp_path): ...
def test_dataset_load_tool_accepts_directory(tmp_path): ...   # 戻り値に adapter と .arf のファイル名
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_dataset_opening.py -q`
Expected: FAIL

- [ ] **Step 3: 実装する**（上の Interfaces のとおり。`handle_tools._load` と `_signature` は `opening` に移し、`handle_tools` からは `opening.open_dataset` を呼ぶ）

- [ ] **Step 4: 通ることを確かめる**（既存のハンドル・読み込み・pipeline の引き継ぎのテストも）

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_dataset_opening.py tests/test_dataset_handle_tools.py tests/test_mztab_tools.py tests/test_dataset_loading_service.py tests/test_readme_links.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add metabolomix/dataset/opening.py metabolomix/mztab/loading.py metabolomix/tools/handle_tools.py metabolomix/tools/mztab_tools.py USAGE.md tests/test_dataset_opening.py
git commit -m "feat(dataset): 解析フォルダ・mzTab-M・job・pipeline を 1 つの入口で開き、job と pipeline の既定を解析フォルダにする"
```

---

### Task 9: 対話のセッションとハンドルの一本化

**Files:**
- Modify: `metabolomix/core/session_state.py`（`AnalysisSession.dataset` / `library` を property に）、`metabolomix/core/handle_dispatch.py`（省略時のロック）、`metabolomix/tools/handle_tools.py`（`dataset_handle_list` に `current`）、`metabolomix/tools/mztab_tools.py`（戻り値にハンドル ID）
- Modify: `tests/conftest.py`（テストごとに新しい登録簿）
- Test: `tests/test_session_handle_unification.py`

**Interfaces:**
- Consumes: `handles.registry`（module 修飾で読む）、`handles.DatasetHandle`、`mztab.polarity.dataset_polarity`
- Produces:
  - `AnalysisSession.current_handle_id: str | None`
  - `AnalysisSession.dataset`（property）: 読むと `registry.peek(current_handle_id)` の `dataset`（無ければ None。閉じられていたら `current_handle_id` も None に戻す）。代入すると、`None` なら `current_handle_id=None`、`DatasetState` なら新しい `DatasetHandle`（`source={"kind": "session"}`、極性は `dataset_polarity`）を `registry.add` して `current_handle_id` にする。
  - `AnalysisSession.library`（property）: 直近のハンドルがあればその `library`、無ければセッションの `LibraryState`。代入も同じ先へ。
  - `handles` と `session_state` の循環は、property の中の遅延 import で避ける。
  - `handle_dispatch._wrap`: `dataset_handle is None` で `session.current_handle_id` のハンドルがあれば、そのロックの下でワーカースレッドに回す。
  - `dataset_handle_list` の各行に `current: bool`。
  - conftest の autouse fixture `_fresh_handle_registry`: `monkeypatch.setattr(handles, "registry", handles.HandleRegistry())`

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_dataset_load_registers_handle_and_sets_current(tmp_path):
    out = dataset_load(directory=str(folder))
    hid = session_state.session.current_handle_id
    assert hid and hid in out
    assert handles.registry.peek(hid).dataset is session_state.session.dataset

def test_closing_current_handle_does_not_fall_back(tmp_path):        # Review Focus 3
    dataset_load(directory=str(folder_a)); first = session_state.session.current_handle_id
    other = <dataset_handle_open(directory=str(folder_b)) のハンドル ID>
    dataset_handle_close(first)
    assert session_state.session.dataset is None
    assert '"missing_state"' in dataset_status()
    assert handles.registry.peek(other) is not None

def test_library_follows_current_handle(): ...
def test_handle_limit_refuses_instead_of_evicting(monkeypatch): ...  # 上限 2 で 3 回目の dataset_load → DATASET_HANDLE_LIMIT、既存の 2 つは開いたまま
def test_interactive_call_waits_for_agent_call_on_same_handle(): ...
    # 既存の handle_dispatch のテスト（`tests/test_dataset_handles.py`）と同じ作りで、
    # 省略呼び出しと明示呼び出しが同じハンドルで重ならないこと
def test_handle_list_marks_current(): ...
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_session_handle_unification.py -q`
Expected: FAIL

- [ ] **Step 3: 実装する**（上の Interfaces のとおり）

- [ ] **Step 4: 全テストを通す**（`session.dataset = ds` を使う既存のテストが多いので全数）

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`
Expected: PASS（落ちたテストは、property の挙動に合わせてテストを直すのではなく、実装が今の意味を保てているかを先に疑う）

- [ ] **Step 5: CLAUDE.md の「データセットのハンドル」の段落を直す**（対話のセッションは直近のハンドルを指すだけで、省略の呼び出しも直近のハンドルのロックの下で動く）

- [ ] **Step 6: コミット**

```bash
git add metabolomix/core/session_state.py metabolomix/core/handle_dispatch.py metabolomix/tools/handle_tools.py metabolomix/tools/mztab_tools.py tests/conftest.py tests/test_session_handle_unification.py CLAUDE.md
git commit -m "refactor(core): 対話のセッションのデータセットをハンドルの登録簿に一本化する"
```

---

### Task 10: 下流は覚えたライブラリを使う

**Files:**
- Modify: `metabolomix/tools/handle_tools.py`（`_open_library` を記憶から開くように）、`metabolomix/library/tools.py`（`library_load` の引数なし・`not_evaluable`・LBM の印）
- Modify: `docs/output_format/library.md`、`USAGE.md`（`library_load` と `dataset_handle_open` の `library` の意味）
- Test: `tests/test_library_from_record.py`

**Interfaces:**
- Consumes: Task 6 `verify_record` / `public_view` / `is_lbm_annotator`、`library.store.open_store(path, read_only=True)`
- Produces:
  - `handle_tools._open_library_for(ds) -> tuple[store | None, str | None]`（記憶が `known` で `verify_record` が通れば開く。`mismatch` / `unknown` / `lbm_only` は開かず理由を返す。極性の設定には落ちない）。`dataset_handle_open(library="auto")` はこれを使う。
  - `library_load()`（引数なし）: 直近のデータセット（または `dataset_handle`）に記憶があればそれを開く。記憶が無ければ今の `resolve_library_path` の順。`file_path` を明示したときは記憶と違っても開き、戻り値に「記憶と違うライブラリ」を書く。
  - `library_match_feature` / `library_msms_qualify`: ライブラリが開かれておらず記憶が `lbm_only` なら `{"status": "not_evaluable", "reason": ...}` を返す（`missing_state` ではない）。記憶が `lbm_also` で、照合した特徴の `annotator_id` に `is_lbm_annotator` が True なら、結果に `identified_by_lbm: true` を足す（None なら足さない）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_handle_opens_remembered_library(tmp_path): ...        # *_Loaded.msp2.dbs が開かれ、戻り値はファイル名だけ
def test_library_sha_mismatch_not_opened(tmp_path):            # Review Focus 4
    ...  # 開いた後で dbs を書き換えて別のハンドルを開く
    assert "不一致" in out and handle.library.store is None
def test_no_record_does_not_use_polarity_setting(tmp_path, monkeypatch): ...
def test_library_load_without_args_uses_record(tmp_path): ...
def test_explicit_file_path_overrides_and_is_reported(tmp_path): ...
def test_lbm_only_dataset_returns_not_evaluable(tmp_path): ...
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests/test_library_from_record.py -q`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: 通ることを確かめる**（ライブラリと再解析のテストも）

Run: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`
Expected: PASS。再解析のテストで `library="auto"` が極性の設定を前提にしているものは、fixture のフォルダに `*_Loaded.msp2.dbs` を置くか job の実効メソッドに MSP を宣言して記憶から開く形に直す（テストの主張は変えない）。

- [ ] **Step 5: コミット**

```bash
git add metabolomix/tools/handle_tools.py metabolomix/library/tools.py docs/output_format/library.md USAGE.md tests/test_library_from_record.py
git commit -m "feat(library): 下流の照合はデータセットが覚えたライブラリを既定にする"
```

---

### Task 11: 実データでの照合と記録

**Files:**
- 追跡対象の変更なし（記録は main ツリーの `docs/HISTRY.md` と `docs/task.md`）

- [ ] **Step 1: reanalysis の run を両方のアダプタで開いて照合する**

`C:\Users\yuu18\datasets\reanalysis\` の下の pipeline の run（`pipeline-run.json` があるもの）ごとに、scratchpad のスクリプトで `open_dataset(pipeline_path=..., source="folder")` と `open_dataset(pipeline_path=..., source="mztab")` を開き `compare_states` を回す。run ごとの `ok`・各検査の数・読み込みの所要秒（解析フォルダのアダプタ）を表にする。スクリプトはリポジトリに置かない（`check.py` も使わない）。

- [ ] **Step 2: 食い違いを分類する**

`ok` が偽の run は、差が「食い違い 1」の既知の差（内部標準の正規化）か、アダプタの不具合かを 1 件ずつ確かめる。不具合なら該当の Task のテストに再現を足して直す（その修正は別のコミット）。

- [ ] **Step 3: 記録する**

main ツリーの `docs/HISTRY.md` に、日付の見出しで「B-1 実データ照合」として run ごとの表・既知の差・直した不具合・LBM の印の結論（Task 6 Step 3）を追記する。`docs/task.md` に B-2・B-3 の未着手を残す。研究室のライブラリのパスは書かない。

---

## この計画が扱わないこと（後続の計画）

- B-2: 公開面の作り直し（`load_dataset` への改名と旧 `load_dataset` / `arf_*` / `arf2_*` / `pai2_*` / `eic_*` / `dcl_parser` の削除、`session.arf` 等のスロットの削除、ツールの対応表 §8、文献と知識の削除、設定と環境変数 §7.4、文書と腐敗防止テスト §10）
- B-3: 定量の工程（§6。msdial-workflow の規則の切り出しと GitHub への配置が先に要る）
- C・A・D はそれぞれの設計書から計画を書く。
