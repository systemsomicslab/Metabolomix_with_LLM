# 再解析の不足を埋める（ハンドル・照合と統計の拡張・証拠の束） Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 再解析のエージェントが、対話と同じ MCP のツールをデータセットのハンドルと宣言的な spec で並行に呼び、判定の数字（行列の作り方・検定・特徴の選び方・MS/MS の合否・証拠の束）を版のあるテスト済みの実装から得られるようにする。

**Architecture:** 1 つのサーバのプロセスにハンドルの登録簿（`core/handles.py`）を置き、ハンドルごとに DatasetState の複製と極性のライブラリを持つ。`dataset_handle` を渡した呼び出しだけをワーカースレッドでハンドルのロックの下に回し、渡さない呼び出しは今までどおりイベントループ上で同期実行する（挙動を一切変えない）。数値の新機能は MCP にもセッションにも依存しない純ロジック（`analysis/recipe_matrix.py` ほか）に置き、ツール層は spec を正規化して呼ぶだけにする。

**Tech Stack:** Python 3.14、FastMCP（mcp 1.30）、numpy / scipy / scikit-learn、**statsmodels（新規。LOWESS）**、anyio（スレッドへの退避）、pytest。

**Spec:** `docs/superpowers/specs/2026-10-22-reanalysis-gaps-design.md`（対: reanalysis-study `docs/specs/2026-10-22-team-efficiency-and-tooled-statistics-design.md`）

## Global Constraints

- 既存の MCP のツールの振る舞い・戻り値・既定値・説明の意味を変えない。足すのは opt-in の引数か新しいツールだけ。
- `dataset_handle` を渡さない呼び出しは変更前と同じ結果・同じ実行方式（イベントループ上の同期実行）になる。
- `dataset_handle` を渡した呼び出しは、ハンドルを渡さない呼び出しが使う `session_state.session` の状態を読みも書きもしない。
- CLI の別の入口は作らない。massbank-context と `analysis/export_contract.py` は変えない。
- 研究室の参照ライブラリはパスを戻り値・結果の JSON・追跡対象の文書・テストに書かない（ファイル名と sha256 だけ）。ハンドルのライブラリは極性で `[library] msp_positive` / `msp_negative`（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）だけから引く（`user_config.get_setting`）。
- spec は正規化して既定値を埋めて返し、未知のキーは拒否する（黙って無視しない）。`rule_hash` は正規化した spec からパスの欄を除いて `core.atomic_io.canonical_hash` で作る。入力ファイル（行列・群・共変量・注入順の TSV）は `input_fingerprint` に sha256 で入れる。
- 同じ spec・同じ入力なら同じ JSON（時刻・乱数の ID を結果に入れない。並びは決定的）。
- 戻り値の作法（CLAUDE.md）: `json_payload`、行の並びは TSV（列名 1 回）、float は丸める、全ツールに `structured_output=False` と `ToolAnnotations`、大きな結果は `output_path` へ JSON を書いて要約だけ返す。
- 新しい環境変数は `LIPIDMIX_` で始める（`LIPIDMIX_MAX_DATASET_HANDLES`）。
- 依存の追加は `pyproject.toml` → `uv lock` → `pyproject.toml` 冒頭の `uv export` → コミット（`tests/test_dependency_lock.py`）。
- テストの fixture はテスト自身が作る。テストは worktree で `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。コミットは pre-commit が全テストを回すのでバックグラウンドで実行し出力をファイルへ落とす。

## Review Focus

1. **閉じた・存在しないハンドル**を渡された → `DATASET_HANDLE_NOT_FOUND` の封筒（開いているハンドルの一覧つき）。グローバルなセッションへ黙って落ちない。テスト: Task 2。
2. **同じハンドルへの並行の呼び出し**（2 人の実行役が同じハンドルを使う）→ ハンドルのロックで直列化され、片方の前処理がもう片方の途中に混ざらない。テスト: Task 1。
3. **TSV の試料名がデータセットと合わない**（群・共変量・注入順・外部行列で、欠け・余り・重複）→ 名前を挙げたエラー。黙って落とさない。テスト: Task 4・5。
4. **検定できない特徴**（全欠測・分散 0・前処理で n<2・OLS の階数不足）→ `status=not_testable` と `reason`、BH の母集団から外す。落ちない。テスト: Task 5。
5. **ハンドルの極性のライブラリが設定に無い／極性が決まらない** → ハンドルはライブラリなしで開き、理由を返す。ライブラリのツールは `missing_state`。逆の極性のライブラリを読まない。テスト: Task 1・3。

---

## ファイルの地図

| ファイル | 役割 |
|---|---|
| `metabolomix/core/handles.py`（新） | ハンドルの登録簿。`DatasetHandle`（`dataset` / `library` / `lock` ほか）、上限、解放。session_state 以外の tools を import しない |
| `metabolomix/core/handle_dispatch.py`（新） | 登録済みツールのうち `dataset_handle` を持つものを包み、ハンドルつきの呼び出しだけをスレッドへ退避する |
| `metabolomix/tools/handle_tools.py`（新） | `dataset_handle_open` / `dataset_handle_list` / `dataset_handle_close` |
| `metabolomix/tools/dataset_analysis_tools.py`・`mztab_tools.py`・`library/tools.py`（改） | 既存ツールへ `dataset_handle`（と `detail`・`with_group_labels`）を足す |
| `metabolomix/library/matching.py`（新） | `library_match_feature` の採点の芯を関数に出す（照合・合否・証拠の束が共有） |
| `metabolomix/library/store.py`（改） | 読み取り専用で開く口と `LibraryStore.search` |
| `metabolomix/msdial/alignment_spot.py`（改） | `library/tools._spot_query` を `spot_msms_query` として移す（dcl ツールと共有） |
| `metabolomix/dcl/tools.py`（改） | `dcl_spot_spectrum` / `dcl_diagnostic_ions` |
| `metabolomix/analysis/recipe_matrix.py`（新） | 行列の recipe v2（順序・検出率・KNN・PQN・LOESS・RSD・半最小値） |
| `metabolomix/analysis/spec_statistics.py`（新） | 検定（student / welch / mannwhitney / ols）・FC・BH・共変量・集約（和・比） |
| `metabolomix/analysis/analysis_spec.py`（新） | analysis-spec.v1 の正規化・`rule_hash`・入力の読み込み（データセット／外部 TSV）と実行 |
| `metabolomix/analysis/msms_rules.py`（新） | MS/MS の合否（msms-rules.v1）・診断イオン・中性損失 |
| `metabolomix/analysis/homolog.py`（新） | 同族列の RT の単調性 |
| `metabolomix/analysis/claim_evidence.py`（新） | 証拠の束（claim-evidence.v1） |
| `metabolomix/analysis/mass_accuracy.py`（新） | m/z 帯ごとの質量精度と新しい知見の候補 |
| `metabolomix/tools/reanalysis_tools.py`（新） | `dataset_analysis_run` / `dataset_claim_evidence` / `dataset_mass_accuracy` / `dataset_novel_candidates` / `dataset_homolog_check` / `dataset_feature_relations` / `library_msms_qualify` |
| `metabolomix/pipeline/inputs.py`・`console/execution.py`・`console/method_file.py`（改） | §4 の 9 (a)(b) |
| `server.py`（改） | 新モジュールの import と `handle_dispatch.install(mcp)` |
| docs（改） | `docs/output_format/{mztab,library,dcl}.md`＋新トピック `reanalysis.md`、`docs/workflow/` の新ツール、`USAGE.md`、`CLAUDE.md` 冒頭の規模表記 |

新しいツールは 13（77 → 90）。

---

### Task 1: ハンドルの登録簿・スレッドへの退避・開閉のツール

**Files:**
- Create: `metabolomix/core/handles.py`, `metabolomix/core/handle_dispatch.py`, `metabolomix/tools/handle_tools.py`
- Modify: `metabolomix/library/store.py`（`open_store(..., read_only: bool = False)`）, `metabolomix/tools/dataset_analysis_tools.py`（`_dataset_polarity` を `metabolomix/mztab/polarity.py` の `dataset_polarity(ds)` へ移し、元の名前は別名で残す）, `server.py`
- Test: `tests/test_dataset_handles.py`

**Interfaces:**
- Produces:
  - `handles.DatasetHandle`: 属性 `handle_id: str`（`dh_` + 12 桁 hex）、`dataset`（DatasetState）、`library`（`session_state.LibraryState`）、`lock`（`threading.RLock`）、`polarity: str | None`、`polarity_source: str`、`arf2_path: str | None`、`dcl_path: str | None`、`source: dict`（読み込み引数）、`created_seq: int`、`last_used_seq: int`。`dataset` と `library` の名前は `AnalysisSession` と同じにする（ツールは `holder.dataset` / `holder.library` を読むだけで両方に効く）。
  - `handles.registry`（module 変数。テストは `handles.registry = handles.HandleRegistry()` で差し替える）。`HandleRegistry.add(handle) -> None`（上限超過で `DomainError("DATASET_HANDLE_LIMIT", ..., {"open_handles": [...]})`）、`get(handle_id) -> DatasetHandle`（無ければ `DomainError("DATASET_HANDLE_NOT_FOUND", ..., {"dataset_handle", "open_handles"})`）、`close(handle_id) -> None`（ライブラリの接続を閉じる）、`list() -> list[DatasetHandle]`、`cached_base(key) / put_base(key, ds, stat_signature)`（読み込み済みの基の LRU、最大 4）。
  - `handles.holder_for(dataset_handle: str | None)`: `None` なら `session_state.session`、そうでなければ `registry.get(...)`。
  - `handles.MAX_HANDLES_ENV = "LIPIDMIX_MAX_DATASET_HANDLES"`（既定 8）。
  - `handle_dispatch.install(mcp) -> list[str]`（包んだツール名）。`handle_dispatch.ALWAYS_THREADED = {"dataset_handle_open"}`。
  - ツール `dataset_handle_open(mztab_path=None, job_path=None, pipeline_path=None, allow_incomplete=False, library="auto", arf2_path=None) -> str`、`dataset_handle_list() -> str`、`dataset_handle_close(dataset_handle: str) -> str`。

- [ ] **Step 1: 失敗するテストを書く**（fixture: テスト内で小さな mzTab-M を書き、陽イオンと陰イオンの 2 本を作る。既存の `tests/test_mztab_tools.py` の fixture 生成を流用してよい）

```python
def test_open_two_handles_do_not_touch_global_session(tmp_path, pos_mztab, neg_mztab):
    session_state.session = AnalysisSession()
    a = json.loads(handle_tools.dataset_handle_open(mztab_path=str(pos_mztab), library="none"))
    b = json.loads(handle_tools.dataset_handle_open(mztab_path=str(neg_mztab), library="none"))
    assert a["dataset_handle"] != b["dataset_handle"]
    assert a["polarity"] == "positive" and b["polarity"] == "negative"
    assert session_state.session.dataset is None

def test_same_input_reuses_base_but_not_state(pos_mztab):
    a = open_(pos_mztab); b = open_(pos_mztab)
    assert b["reused_base"] is True
    ha, hb = handles.registry.get(a["dataset_handle"]), handles.registry.get(b["dataset_handle"])
    assert ha.dataset is not hb.dataset          # 複製（片方の前処理がもう片方に見えない）

def test_handle_limit(monkeypatch, pos_mztab):
    monkeypatch.setenv("LIPIDMIX_MAX_DATASET_HANDLES", "1")
    open_(pos_mztab)
    err = json.loads(handle_tools.dataset_handle_open(mztab_path=str(pos_mztab), library="none"))
    assert err["error"]["code"] == "DATASET_HANDLE_LIMIT"

def test_close_then_list(pos_mztab): ...  # close 後は list に無く、get は DATASET_HANDLE_NOT_FOUND

def test_library_auto_without_config_opens_without_library(pos_mztab):
    out = json.loads(handle_tools.dataset_handle_open(mztab_path=str(pos_mztab)))
    assert out["library"] is None and "msp_positive" in out["library_note"]

def test_library_auto_uses_polarity_setting_only(monkeypatch, tmp_path, pos_mztab):
    # 自作の小さな .msp を MSDIAL_MSP_POS に置き、MSDIAL_MSP_NEG は別ファイル
    out = open_(pos_mztab, library="auto")
    assert out["library"]["file"] == "pos.msp" and len(out["library"]["sha256"]) == 64
    assert str(tmp_path) not in json.dumps(out)  # パスを出さない

def test_handle_calls_run_off_loop_and_in_parallel(two_handles):
    # mcp.call_tool 経由。tool A の実体が threading.Event を待ち、tool B がそれを set する。
    # 直列なら A がタイムアウトする。両方が非メインスレッドで走ったことも確かめる。

def test_same_handle_calls_are_serialized(one_handle): ...  # 2 つの呼び出しの区間が重ならない

def test_call_without_handle_stays_on_event_loop_thread(): ...  # 実行スレッド = ループのスレッド
```

- [ ] **Step 2:** `pytest tests/test_dataset_handles.py -q` → 失敗（モジュールが無い）。
- [ ] **Step 3: 実装する**
  - `dataset_handle_open` は `dataset_load` と同じ読み込み関数（`load_dataset_state` / `session_handoff.load_run_dataset`）を呼び、`session_state.session` に代入しない。基のキャッシュのキーは `canonical_hash({"kind", "path": 絶対パス, "allow_incomplete"})`、再利用はキャッシュ時の `ds.source_files` の各パスの `(size, mtime_ns)` が今と同じときだけ。ハンドルには `copy.deepcopy(base)` を渡す。
  - 極性は `dataset_polarity(ds)`。`library="auto"` は極性が決まれば `user_config.get_setting("library.msp_positive"|"library.msp_negative")` だけを見て `open_store(path, read_only=True)`（未構築なら構築してから。構築はプロセス内の 1 本の `threading.Lock` の下）。極性不明・未設定・`library="none"` はライブラリなしで開き `library_note` に理由。`read_only=True` の接続は `sqlite3.connect(f"file:{path}?mode=ro", uri=True, check_same_thread=False)`。
  - `arf2_path` 省略時は `ds.feature_qc["join"]["path"]`（検出状態を取り込んだ `.arf`）から `alignment_spot.alignment_sibling(arf, ".arf2")`、無ければ `curation_apply.arf2_for_mztab`。`.dcl` は `.arf2` の兄弟。
  - 戻り値（json_payload）: `dataset_handle, polarity, polarity_source, source_format, n_features, n_samples, detection_available, library ({file, sha256, record_count} | null), library_note, arf2 (ファイル名), dcl (ファイル名), reused_base, open_handles, max_handles`。
  - `handle_dispatch.install`: `mcp._tool_manager._tools` のうち関数の引数に `dataset_handle` を持つものと `ALWAYS_THREADED` を、`async def runner(**kw)` に包み `tool.is_async = True`。`kw.get("dataset_handle") is None` かつ ALWAYS_THREADED でなければ元の関数をそのまま同期で呼ぶ。そうでなければ `anyio.to_thread.run_sync` の中で `registry.get(h).lock` を取って呼ぶ（未知のハンドルはロックを取らずにそのまま呼び、ツール側の封筒に任せる）。`server.py` の末尾で 1 回呼ぶ。
- [ ] **Step 4:** テストが通る。`tests/test_server_registration.py` の EXPECTED_TOOLS に 3 本を足して通す（件数のコメントも 77→80）。
- [ ] **Step 5:** コミット `feat(handles): データセットのハンドルを開閉し、ハンドルつきの呼び出しをスレッドで回す`

### Task 2: 既存のツールへ `dataset_handle` を足す

**Files:**
- Modify: `metabolomix/tools/dataset_analysis_tools.py`（`dataset_preprocess` `dataset_pca` `dataset_differential` `dataset_export_differential` `dataset_export_evidence_bundle` `dataset_set_sample_metadata` `dataset_statistic` `dataset_build_matrix` `dataset_feature_search`）、`metabolomix/tools/mztab_tools.py`（`dataset_status`）、`metabolomix/library/tools.py`（`library_load` `library_match_feature` `library_plot_mirror`）
- Test: `tests/test_dataset_handle_tools.py`

**Interfaces:**
- Consumes: `handles.holder_for`, `handles.registry`（Task 1）
- Produces: 上の 13 ツールの末尾に `dataset_handle: str | None = None`。共通ヘルパ `metabolomix/tools/_holder.py: resolve_holder(dataset_handle) -> holder | str`（str は `mztab_error("DATASET_HANDLE_NOT_FOUND", ...)` の封筒）。ツールの本文は `session_state.session.dataset` / `.library` を `holder.dataset` / `holder.library` に置き換えるだけ（`dataset_handle=None` では `holder is session_state.session`）。

既存の `dcl_parser` / `dcl_find_msms` はファイルのパスで動きセッションのデータセットを読まないので、`dataset_handle` を足さない（ハンドルから `.dcl` を引く口は Task 3・7 の新しい dcl ツールが持つ）。

- [ ] **Step 1: 失敗するテスト**

```python
@pytest.mark.parametrize("tool,kwargs", ALL_HANDLE_TOOLS)   # 13 本と最小の引数
def test_unknown_handle_is_an_error_not_global(tool, kwargs):
    session_state.session.dataset = make_ds()               # グローバルにはデータがある
    out = json.loads(tool(**kwargs, dataset_handle="dh_000000000000"))
    assert out["error"]["code"] == "DATASET_HANDLE_NOT_FOUND"

def test_preprocess_on_handle_leaves_global_untouched(handle, global_ds):
    dataset_preprocess(normalize="pqn", dataset_handle=handle)
    assert global_ds.pp_matrix is None
    assert handles.registry.get(handle).dataset.pp_matrix is not None

def test_without_handle_output_is_byte_identical(global_ds):
    # 変更前の戻り値を固定した文字列と比べる（dataset_status / dataset_feature_search の窓検索 /
    # dataset_differential の要約）。server_version など時刻依存の欄は事前に固定する。

def test_library_match_with_handle_uses_handle_store_and_arf2(handle_with_library, fixture_arf2_dcl): ...
def test_plot_mirror_reads_handle_last_match(handle_with_library): ...   # グローバルの last_match は None のまま
```

- [ ] **Step 2:** 失敗を確かめる。
- [ ] **Step 3:** 置き換える。`library_match_feature` は `spot_id` の解決で `arf2_path is None and dataset_handle` のとき `holder.arf2_path` と `holder.dcl_path` を使う。docstring に「`dataset_handle`: `dataset_handle_open` が返したハンドル。省略時は対話のセッション。再解析の並行のエージェント向け」を 1 行足す（既存の文は変えない）。
- [ ] **Step 4:** 全テスト（既存の `tests/test_dataset_*` `tests/test_library_*` が無変更で通ること）。
- [ ] **Step 5:** コミット `feat(handles): dataset_* と library_* がハンドルを受ける（省略時は従来どおり）`

### Task 3: 照合の戻り値の拡張・ライブラリの検索・spot_id の代表スペクトル

**Files:**
- Create: `metabolomix/library/matching.py`
- Modify: `metabolomix/library/tools.py`, `metabolomix/library/store.py`, `metabolomix/analysis/spectral_match.py`, `metabolomix/msdial/alignment_spot.py`, `metabolomix/dcl/tools.py`, `docs/output_format/library.md`（`matched_peaks_count` の前駆体の扱い）
- Test: `tests/test_library_detail.py`

**Interfaces:**
- Produces:
  - `spectral_match.fragment_matches_excluding_precursor(alignment: list[dict], precursor_mz: float, tol: float) -> int | None`（`alignment` は `match_spectrum` の `alignment`。`matched` かつ `abs(mz - precursor_mz) > tol` の数。alignment が空なら None）。
  - `matching.score_candidates(store, *, precursor_mz, rt, ion_mode, measured, mz_tol, ms2_tol, rt_tol) -> dict`：`{"status": "success"|"no_candidates", "query": {...}, "scoring": {...}, "scored": [(record, result), ...]（total_score 降順・安定）, "unscored_mz", "n_candidates"}`。`library_match_feature` はこれを呼ぶ形に直す（戻り値は 1 文字も変えない）。
  - `alignment_spot.spot_msms_query(spot_id, *, arf2_path, dcl_file, rt, ion_mode) -> dict`（今の `library/tools._spot_query` をそのまま移し、元の名前は別名で残す）。
  - `library_match_feature(..., detail: bool = False, precursor_exclusion_tol: float | None = None)`：`detail=True` のときだけ TSV の末尾に `inchikey` `fragment_matches_excl_precursor` 列、payload に `precursor_exclusion_tol`（省略時は解決済みの `ms2_tol`）、`measured_peaks`（`[[mz(4桁), 強度(有効4桁)]]`）、`alignments`（`[{"rank", "matched": [[mz, measured, reference]]}]`、一致した行だけ）。
  - `LibraryStore.search(*, name=None, formula=None, inchikey_block=None, ion_mode=None, limit=50) -> tuple[list[dict], int]`（name は大小無視の部分一致、formula は完全一致、inchikey は先頭 14 文字の一致、ion_mode は NOCASE。返す件数と打ち切った件数）。
  - ツール `library_search(name=None, formula=None, inchikey=None, ion_mode=None, limit=50, dataset_handle=None)`：TSV 列 `name precursor_mz ion_mode adduct rt formula inchikey ontology compound_class n_peaks`。条件が 1 つも無ければ `invalid_arguments`。
  - ツール `dcl_spot_spectrum(spot_id: int, arf2_path=None, dcl_file=None, dataset_handle=None)`：`{status, spot_id, precursor_mz, rt, ion_mode, n_peaks, peaks}`（失敗は `spot_msms_query` の封筒そのまま）。

- [ ] **Step 1: 失敗するテスト**

```python
def test_fragment_count_excludes_precursor():
    alignment = [{"mz": 184.07, "measured": 1, "reference": 1, "matched": True},
                 {"mz": 760.585, "measured": 1, "reference": 1, "matched": True}]
    assert fragment_matches_excluding_precursor(alignment, 760.5851, 0.01) == 1

def test_detail_false_output_unchanged(lib_store, fixture_spot):   # 既存の戻り値と完全一致
def test_detail_true_adds_columns_and_alignments(lib_store, fixture_spot):
    out = json.loads(library_match_feature(spot_id=1, detail=True, ...))
    assert out["candidates_table"].splitlines()[0].endswith("inchikey\tfragment_matches_excl_precursor")
    assert out["alignments"][0]["rank"] == 1 and "measured_peaks" in out
def test_search_by_name_formula_inchikey(lib_store): ...
def test_search_store_is_read_only(tmp_path, small_msp):   # read_only=True の接続で INSERT が失敗する
def test_spot_spectrum_returns_peaks(fixture_arf2_dcl): ...
```

- [ ] **Step 2–4:** 失敗 → 実装 → 通る。`docs/output_format/library.md` の `matched_peaks_count` に「前駆体の一致も 1 本に数える。前駆体を除いた本数は `detail=true` の `fragment_matches_excl_precursor`」を足し、新しい列・ツールの欄を書く。
- [ ] **Step 5:** コミット `feat(library): 照合の詳細（前駆体を除いた断片の本数）・ライブラリの検索・spot の代表スペクトル`

### Task 4: 行列の recipe v2

**Files:**
- Create: `metabolomix/analysis/recipe_matrix.py`
- Modify: `pyproject.toml`（`statsmodels>=0.14`）, `uv.lock`, `requirements*.lock.txt`
- Test: `tests/test_recipe_matrix.py`

**Interfaces:**
- Produces:
  - `RECIPE_SCHEMA = "matrix-recipe.v2"`、`normalize_recipe(recipe: dict) -> dict`（既定を埋める・未知のキーと「`steps` に無いのに設定がある段」を `DomainError("MATRIX_RECIPE_INVALID")` で拒否）。
  - `build_matrix(values: np.ndarray, sample_ids: list[str], feature_ids: list[str], roles: dict[str, str], recipe: dict, *, detected: np.ndarray | None = None, injection_order: dict[str, float] | None = None) -> dict`：`values` と `detected` は（試料 × 特徴）。戻り値 `{"values", "sample_ids", "feature_ids"（残った特徴）, "steps": [段ごとの報告], "n_features_in", "n_features_out"}`。入力の配列は書き換えない。
  - recipe の形（正規化後。既定値はこのとおり）:

```json
{"schema": "matrix-recipe.v2",
 "steps": [],
 "detection_filter": {"basis": "detected", "min_rate": {}},
 "impute": {"method": "half_min", "missing": "zero", "k": 10, "scale": "log2", "rowmax": 0.5, "colmax": 0.8},
 "normalize": {"method": "pqn", "reference": "qc"},
 "drift_correct": {"method": "loess", "span": 0.75, "iterations": 0, "min_qc": 4},
 "qc_rsd_filter": {"max_rsd": 0.3}}
```
  `steps` は `detection_filter` `impute` `normalize` `drift_correct` `qc_rsd_filter` の並び（重複不可）。設定の塊は `steps` に挙げた段の分だけ正規化後に残す（挙げない段の塊が入力にあれば拒否）。`steps: []` は「前処理なし」。

  定義（docs に写す）:
  - 検出の基準 `basis`: `detected` は gap-fill を除く実検出（`detected` が None なら `DomainError("DETECTION_UNAVAILABLE")`）、`positive` は有限かつ > 0。`min_rate` は役割（`qc` `sample` `blank`）ごとの最小の率で、挙げた役割のすべてを満たす特徴を残す。
  - 欠測 `missing`: `zero` は ≤0 か非有限、`nan` は非有限だけ。補完の後、欠測でなかったセルは変わらない。
  - `half_min`: 特徴ごとの正の有限値の最小の半分（正の値が無い特徴は欠測のまま）。`none`: 補完しない。
  - `knn`: R `impute::impute.knn`（`maxp` = 特徴数、すなわち 2 段の分割をしない）と同じ定義。特徴を行に、`scale="log2"` なら正の値の log2 の上で、欠測のある特徴ごとに、その特徴が観測された列のうち相手も観測された列での平方差の平均を距離にして（共通の列が無い相手は除く）近い `k` 個を選び（同距離は特徴の並び順）、欠測の列をその近傍の観測値の平均で埋める。近傍がすべて欠測ならその列（試料）の平均。欠測の割合が `rowmax` を超える特徴は列の平均で埋める。欠測の割合が `colmax` を超える試料があれば `DomainError("RECIPE_PRECONDITION")`。距離はブロック化した行列積で求める（33,000 特徴 × 60 試料で数秒）。
  - `normalize`: `pqn` は参照（`reference` の役割の試料の、特徴ごとの正の有限値の中央値）に対する商（両方が正の有限な特徴だけ）の中央値で各試料を割る。`tic` は有限値の和、`median` は正の有限値の中央値を試料間の中央値で割った係数。係数が 0・非有限の試料は割らずに報告する。
  - `drift_correct`: 特徴ごとに、正の有限な QC 値を注入順に対して `statsmodels.nonparametric.smoothers_lowess.lowess(endog, exog, frac=span, it=iterations, delta=0.0, xvals=全試料の注入順)` で当て、`値 / 当てはめ × QC 値の中央値` に直す。正の QC 値が `min_qc` 未満、または当てはめが正の有限でない特徴は直さず数を報告する。注入順はすべての試料に要る（欠ければ名前を挙げて `DomainError("INJECTION_ORDER_MISSING")`）。
  - `qc_rsd_filter`: QC の有限値の `std(ddof=1)/mean` が `max_rsd` 以下の特徴を残す（RSD が出ない特徴は落とし、別に数える）。

- [ ] **Step 1: 失敗するテスト**（すべて手で値を計算できる小行列で）

```python
def test_steps_run_in_given_order():  # 同じ段でも順序で値が変わる例（impute→normalize と normalize→impute）
def test_unknown_key_and_unlisted_block_rejected():
def test_defaults_are_filled_and_echoed():
def test_detection_filter_by_role_positive_vs_detected():   # gap-fill のセルが率に入るかで結果が分かれる
def test_knn_matches_hand_computed_impute_knn():            # 5 特徴 × 4 試料、k=2。期待値は手計算（R impute.knn と同じ）
def test_knn_rowmax_and_colmax():
def test_pqn_reference_qc():
def test_loess_matches_statsmodels_direct_call():           # 期待値を lowess(...xvals=...) で直接計算して比べる
def test_drift_requires_order_for_every_sample():
def test_qc_rsd_filter_drops_nonfinite_rsd():
def test_inputs_not_mutated():
```

- [ ] **Step 2:** 失敗を確かめる。
- [ ] **Step 3:** 依存を足す（`uv add statsmodels` 相当の手順）→ 実装。
- [ ] **Step 4:** 通る（`tests/test_dependency_lock.py` も）。
- [ ] **Step 5:** コミット `feat(analysis): 行列の recipe v2（順序・役割ごとの検出率・特徴方向の KNN・LOESS）`

### Task 5: 統計の spec（検定・FC・補正・共変量・集約・感度解析・外部の行列）と `dataset_analysis_run`

**Files:**
- Create: `metabolomix/analysis/spec_statistics.py`, `metabolomix/analysis/analysis_spec.py`, `metabolomix/tools/reanalysis_tools.py`
- Test: `tests/test_spec_statistics.py`, `tests/test_analysis_spec.py`

**Interfaces:**
- Consumes: `recipe_matrix.normalize_recipe` / `build_matrix`（Task 4）、`handles.holder_for`
- Produces:
  - `spec_statistics.test_feature(ref: np.ndarray, test: np.ndarray, *, kind, transform, design=None) -> dict`：`{"status", "reason", "statistic", "p", "beta_log2"}`。kind は `student`（等分散 t）・`welch`・`mannwhitney`（両側、`scipy.stats.mannwhitneyu(..., alternative="two-sided")`、変換に依らない）・`ols`（`transform(y) ~ 1 + group(test=1) + 共変量`、群の係数の t の両側 p。階数不足・自由度 0 は not_testable）。`transform="log2"` は正の有限値だけ（それ以外は欠測）。各群の有限値が 2 未満・両群とも分散 0 は not_testable。
  - `spec_statistics.fold_change(ref_raw, test_raw, definition) -> float | None`：`arithmetic_mean_ratio`（変換前の算術平均の比 test/ref）か `geometric_mean_ratio`（log2 の平均の差の 2 乗）。平均が 0 以下なら None。
  - `spec_statistics.bh(p: list[float | None]) -> list[float | None]`（None は母集団に入れない）。
  - `analysis_spec.ANALYSIS_SCHEMA = "analysis-spec.v1"`、`normalize_spec(spec: dict) -> dict`、`rule_hash(normalized: dict) -> str`、`run_spec(normalized: dict, ds=None) -> dict`（結果 `analysis-result.v1`）。
  - ツール `dataset_analysis_run(spec: dict, dataset_handle: str | None = None, output_path: str | None = None) -> str`。

  spec の形（正規化後。`|` は択一、既定は右の値）:

```json
{"schema": "analysis-spec.v1",
 "input": {"source": "dataset"} | {"source": "tsv", "matrix_path": "...", "orientation": "features_by_samples", "feature_column": "feature_id"},
 "groups": {"source": "metadata"} | {"source": "tsv", "path": "...", "sample_column": "sample", "group_column": "group", "role_column": null},
 "comparison": {"reference": "HC", "test": "RA"},
 "injection_order": null | {"source": "metadata"} | {"source": "tsv", "path": "...", "sample_column": "sample", "order_column": "order"},
 "covariates": null | {"path": "...", "sample_column": "sample", "missing": "drop",
                       "columns": [{"name": "sex", "type": "categorical", "reference": "M"}, {"name": "age", "type": "numeric"}]},
 "features": {"all": true} | {"feature_ids": ["1735"]},
 "aggregate": null | {"kind": "sum", "label": "PC_total", "feature_ids": [...]} | {"kind": "sum", "label": "...", "select": {"name_regex": "^PC ", "evidence": ["sme"]}}
                   | {"kind": "ratio", "label": "...", "numerator": "1735", "denominator": "1244"},
 "analyses": [{"analysis_id": "main", "matrix": <matrix-recipe.v2>, "test": {"kind": "student", "transform": "log2", "covariates": []},
               "fold_change": "arithmetic_mean_ratio", "correction": "bh", "bh_population": "matrix",
               "thresholds": {"p": 0.05, "q": null, "fc": 1.5}}]}
```
  - 群は `metadata`（`ds.sample_metadata_rows` の `group`・`role`・`include`。`include=true` かつ `role=sample` だけを群に入れる。メタデータが無ければ `DomainError("GROUPS_UNAVAILABLE")`）か TSV。外部の行列（`input.source=tsv`）は群も TSV で与える（`role_column` があれば役割、無ければ全員 `sample`）。行列の試料は群の TSV に無い試料も役割 `sample`・群なしとして持つ（QC は役割で区別）。
  - 試料名の照合は完全一致。群・注入順・共変量の TSV に行列の試料が欠ける／余る／重複するときは名前を挙げて `DomainError("SAMPLE_MISMATCH")`（共変量の欠測値は `missing: drop` なら OLS の解析でその試料を除き `n` を報告、`error` なら拒否）。
  - `aggregate` は各解析の行列を作った**後**に、試料ごとに選んだ特徴の値の和（欠測は欠測のまま。和を取る特徴のどれかが欠測の試料は欠測）か比を作り、1 本の特徴（id = `label`）として検定する。行列で落ちた特徴は和から外して `aggregate.dropped` に挙げる。`select` は SME の rank 1 の名前（`evidence` に `sml` があれば SML の名前も）に正規表現を当てる。選んだ特徴の一覧を結果に返す。
  - `bh_population`: `matrix`（その解析で検定できた全特徴）か `reported`（報告する特徴だけ）。`correction: none` なら q は null。
  - 行の欄（解析ごと）: `feature_id status reason group_ref group_test n_ref n_test mean_ref mean_test fc_test_over_ref log2fc_test_over_ref statistic p q beta_log2 direction meets_thresholds`。`direction` は `test_higher` / `ref_higher`（ols は `beta_log2` の符号、ほかは FC）。`meets_thresholds` は与えた閾値（null の閾値は見ない）をすべて満たすか。**再現したかの判定は返さない。**
  - 結果: `{"schema": "analysis-result.v1", "spec", "rule_hash", "input_fingerprint": {"dataset_hash", "files": {role: sha256}}, "comparison": "RA_vs_HC", "group_ref", "group_test", "aggregate", "analyses": [{"analysis_id", "matrix": {"n_features_in", "n_features_out", "steps"}, "bh_population_n", "thresholds", "rows"}]}`。`rule_hash` は `spec` からパスの欄（`matrix_path` `path`）を除いて作る。
  - ツールの戻り値: `output_path` があれば結果の JSON を書き（`ensure_ascii=False`、`sort_keys=True`、改行なし）、要約（`rule_hash`・解析ごとの件数・`output_path`）だけ。無ければ報告する行が解析あたり 50 以下のときだけ行を TSV で返し、超えれば `DomainError("ANALYSIS_OUTPUT_REQUIRED")`。行の数値は p・q・平均を有効 4 桁、FC を 4 桁に丸める（JSON のファイルは丸めない）。

- [ ] **Step 1: 失敗するテスト**

```python
def test_student_log2_and_arithmetic_fc_match_scipy():      # scipy.stats.ttest_ind(log2, equal_var=True) と一致
def test_mannwhitney_two_sided():
def test_ols_group_p_matches_closed_form():                  # numpy の正規方程式で作った期待値
def test_ols_rank_deficient_is_not_testable():
def test_not_testable_excluded_from_bh():
def test_bh_population_matrix_vs_reported():
def test_group_names_and_direction_on_every_row():           # comparison == "RA_vs_HC"、group_ref == "HC"
def test_sensitivity_analyses_return_each_numbers():         # 前処理なし / half_min / welch / mannwhitney / ols の 5 本
def test_same_spec_same_json_and_rule_hash(tmp_path):        # 2 回走らせて output_path の JSON がバイト一致
def test_rule_hash_ignores_paths_but_fingerprint_tracks_files(tmp_path):
def test_external_matrix_tsv_input(tmp_path):                # データセット無しで動く（較正の口）
def test_sample_mismatch_names_listed(tmp_path):
def test_covariate_missing_drop_reports_n(tmp_path):
def test_aggregate_sum_lists_features_and_dropped():
def test_aggregate_ratio():
def test_unknown_key_rejected():
def test_handle_input_does_not_touch_global(handle): ...
```

- [ ] **Step 2–4:** 失敗 → 実装 → 通る。
- [ ] **Step 5:** コミット `feat(analysis): 統計の spec（Student・Mann-Whitney・共変量の OLS・集約・感度解析・外部の行列）`

### Task 6: 群の向きを統計の各行に（`dataset_feature_search` の opt-in）

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`（`STATS_COLUMNS_LABELED`）, `metabolomix/tools/dataset_analysis_tools.py`
- Test: `tests/test_feature_search_labels.py`

**Interfaces:**
- Produces: `dataset_feature_search(..., with_group_labels: bool = False)`。True のときだけ `stats` の TSV の末尾に `comparison group_a group_b log2fc_test_over_ref`（`comparison` は `provenance.comparison.comparison_id`、無ければ `<group_b>_vs_<group_a>`。`log2fc_test_over_ref` は v1 の `log2fc` と同じ値＝ group_b / group_a）。

- [ ] **Step 1:** テスト: False の戻り値が変更前と一致 / True で 4 列が足され、`RA_vs_HC`・`group_a=HC` の例で `log2fc_test_over_ref == log2fc`。
- [ ] **Step 2–4:** 失敗 → 実装 → 通る。`docs/output_format/mztab.md` に列を足す。
- [ ] **Step 5:** コミット `feat(feature_search): 統計の行に群の名前と比較の向きの log2FC（opt-in）`

### Task 7: MS/MS の合否・診断イオン・同族列

**Files:**
- Create: `metabolomix/analysis/msms_rules.py`, `metabolomix/analysis/homolog.py`
- Modify: `metabolomix/tools/reanalysis_tools.py`（`library_msms_qualify` `dataset_homolog_check`）, `metabolomix/dcl/tools.py`（`dcl_diagnostic_ions`）
- Test: `tests/test_msms_rules.py`, `tests/test_homolog.py`

**Interfaces:**
- Consumes: `matching.score_candidates`、`spectral_match.fragment_matches_excluding_precursor`、`alignment_spot.spot_msms_query`（Task 3）
- Produces:
  - `msms_rules.RULES_SCHEMA = "msms-rules.v1"`、`normalize_rules(rules) -> dict`、`qualify(scored: list[tuple[dict, dict]], measured, precursor_mz, rules) -> dict`、`find_ions(measured, ions, *, precursor_mz) -> list[dict]`。
  - rules の形（正規化後）:

```json
{"schema": "msms-rules.v1",
 "target": {"inchikey_block": "HVYWMOMLDIMFJA", "formula": null},
 "top_rank": 5,
 "min_fragments_excl_precursor": 2,
 "precursor_exclusion_tol": null,
 "same_formula_margin": null | {"score": "total_score", "min_delta": 0.0},
 "lipid_class": null | "PC",
 "diagnostic_ions": [{"mz": 184.0733, "label": "phosphocholine", "tol": 0.01, "min_relative_intensity": 0.0}],
 "neutral_losses": [{"delta": 141.0191, "label": "PE head", "tol": 0.01}],
 "ions_match": "all"}
```
  - 検査（`checks` の各行 `{check, status: pass|fail|not_evaluable, value, detail}`）: `target_in_top_rank`（`top_rank` 位までに InChIKey 先頭ブロックが一致する候補がある。最上位のそれを「目標の候補」）、`min_fragments`（目標の候補の前駆体を除いた一致の本数 ≥ 最小）、`same_formula_margin`（全候補のうち目標と同じ分子式で先頭ブロックが違うものの最高点との差 ≥ `min_delta`。別物が無ければ pass）、`lipid_class`（1 位の候補の名前を `normalize_lipid_name` で読んだクラス＝正規化名の最初の空白の前が一致）、`diagnostic_ions`・`neutral_losses`（測定スペクトルに ±tol のピークがあり、相対強度 ≥ 最小。`ions_match` が `all` なら全部、`any` なら 1 つ）。設定していない検査は行を出さない。`verdict`: すべて pass → `qualified`、どれか fail → `not_qualified`、fail が無く not_evaluable がある・MS/MS が無い → `not_evaluable`。
  - ツール `library_msms_qualify(rules: dict, spot_id=None, precursor_mz=None, rt=None, ion_mode=None, arf2_path=None, dcl_file=None, top_n: int = 10, dataset_handle=None) -> str`：`{verdict, checks, rule_hash, rules, candidates (TSV: rank name formula inchikey adduct total_score matched_peaks_count fragment_matches_excl_precursor), query}`。
  - ツール `dcl_diagnostic_ions(spot_ids: list[int], ions: list[dict], arf2_path=None, dcl_file=None, dataset_handle=None) -> str`：ions の各要素は `{"mz"|"neutral_loss", "label", "tol"=0.01}`。spot は 200 まで。TSV `spot_id label kind target_mz found observed_mz relative_intensity status`（status: `ok` / spot の失敗理由）。
  - `homolog.check(members: list[dict], rule: dict) -> dict`、ツール `dataset_homolog_check(spec: dict, dataset_handle=None) -> str`。spec:

```json
{"members": [{"feature_id": "1", "carbons": 34, "double_bonds": 1}] | null,
 "auto": null | {"lipid_class": "PC", "evidence": ["sme"], "adduct": null},
 "rule": {"carbon_direction": "increasing", "double_bond_direction": "decreasing", "max_violations": 0}}
```
  （members と auto はどちらか一方。auto は SME rank 1 の名前を `normalize_lipid_name` で種のレベル `CLASS C:DB` に読めた特徴を集める。）同じ（C, DB）の特徴が複数あれば RT の中央値で並べ、`duplicates` に挙げる。同じ DB の中で C の昇順に並べた隣どうし、同じ C の中で DB の昇順に並べた隣どうしで、規則の向きに反する対を違反として数える。戻り値 `{verdict: pass|fail, n_violations, violations (TSV), members (TSV: feature_id carbons double_bonds rt status), duplicates, rule_hash}`。

- [ ] **Step 1: 失敗するテスト**（スペクトルと候補は手で作る）: 前駆体だけの一致は `min_fragments` で fail / 同じ分子式の別物が上なら margin で fail / 診断イオンの有無 / 中性損失 / 未設定の検査は出ない / MS/MS 無しは not_evaluable / 規則の既定と並びに依らない rule_hash / 同族列の違反の数え方と duplicates / auto が名前から集める。
- [ ] **Step 2–4:** 失敗 → 実装 → 通る。
- [ ] **Step 5:** コミット `feat(analysis): MS/MS の合否・診断イオン・同族列の RT の検査`

### Task 8: 特徴どうしの関係

**Files:**
- Modify: `metabolomix/tools/reanalysis_tools.py`（`dataset_feature_relations`）
- Test: `tests/test_feature_relations.py`

**Interfaces:**
- Consumes: `curation.relations.find_relations`、`profile_correlation`
- Produces: `analysis/claim_evidence.feature_relations(ds, feature_id, *, rt_window, mz_tol, min_r, adducts) -> dict`（`{"coeluting": [...], "relations": [...]}`。相手は `ds.feature_ids` の全特徴で、`spot_id` に SMF_ID、`rt` は分、`heights` は `feature_matrix` の行、`name`・`formula` は SME rank 1。MS-DIAL のリンクは無い（`links={}`）ので裏付けは強度の相関だけ）。ツール `dataset_feature_relations(feature_id: str, rt_window: float = 0.05, mz_tol: float = 0.01, min_r: float = 0.7, adducts: list[str] | None = None, dataset_handle=None)`：`adducts` 省略時はデータセットの極性の `DEFAULT_ADDUCT_TIERS` を平らにしたもの。TSV 2 本（`coeluting`: `feature_id mz rt drt r name`、`relations`: `of relation dmz_mda drt r isotope_consistent strong`）。

- [ ] **Step 1:** テスト: [M+H]+ と [M+Na]+ の組・M+1 同位体・-H2O の源内断片を持つ 4 特徴の fixture で、関係が名前つきで返る／RT 窓の外は出ない／未知の feature_id はエラー。
- [ ] **Step 2–4:** 失敗 → 実装 → 通る。
- [ ] **Step 5:** コミット `feat(analysis): mzTab の特徴どうしの関係（同時溶出・付加体・同位体・源内断片）`

### Task 9: 証拠の束

**Files:**
- Modify: `metabolomix/analysis/claim_evidence.py`, `metabolomix/tools/reanalysis_tools.py`（`dataset_claim_evidence`）
- Test: `tests/test_claim_evidence.py`

**Interfaces:**
- Consumes: `feature_search.normalize_targets` / `match_target` / `rule_hash`、`msms_rules`、`matching.score_candidates`、`feature_relations`
- Produces: `claim_evidence.build(ds, spec, *, store, arf2_path, dcl_path, polarity) -> dict`、ツール `dataset_claim_evidence(spec: dict, dataset_handle=None, output_path=None) -> str`。spec:

```json
{"claim_id": "C1",
 "target": {"target_id": "C1", "inchikey": "...", "formula": "...", "rt_min": 8.4},
 "scan": {"ppm": [10.0, 20.0], "rt_tol_min": [0.3, 1.0], "adduct_tiers": null},
 "msms": null | {"top_n": 10, "rules": {<msms-rules.v1 から target を除いたもの>}},
 "relations": null | {"rt_window": 0.05, "mz_tol": 0.01, "min_r": 0.7},
 "max_candidates": 50}
```
  結果 `claim-evidence.v1`: `{claim_id, rule_hash, input_fingerprint: {dataset_hash, library_sha256}, rule, default_selection: {status, reason, selected_feature_id, selected_adduct, selected_evidence, counts}（target をそのまま match_target した結果＝ dataset_feature_search の目標照合と同じ規則）, scans: [{ppm, rt_tol_min, status, selected_feature_id, candidate_feature_ids}]（ppm × rt_tol_min の全組）, candidates: [{feature_id, mz, rt, evidence, annotated_name, adducts, msms: {status, top (行の list), qualify}, relations}]（全走査の候補の和を feature_id で並べ、max_candidates で打ち切り）, truncated, notes}`。統計は入れない（実行役が選んだ特徴で `dataset_analysis_run` を別に呼ぶ）。MS/MS は spot_id = SMF_ID として `spot_msms_query`（前駆体の食い違いは `msms.status` に残す）。ライブラリが無いとき `msms` は `{"status": "no_library"}`。
  戻り値: `output_path` があれば JSON を書いて要約（`rule_hash` `default_selection` 候補数）、無ければ全体（候補が 10 を超えれば `output_path` を求める）。

- [ ] **Step 1:** テスト: 2 回走らせて JSON がバイト一致・`rule_hash` 一致 / 既定の特徴が `dataset_feature_search` の目標照合と同じ / 走査を広げると候補が増える / ライブラリ無しでも束は作れる / 解決できないと `needs_input`。
- [ ] **Step 2–4:** 失敗 → 実装 → 通る。
- [ ] **Step 5:** コミット `feat(analysis): 1 主張の証拠の束（候補の走査・照合・合否・関係・既定の特徴）`

### Task 10: 質量精度と新しい知見の候補

**Files:**
- Create: `metabolomix/analysis/mass_accuracy.py`
- Modify: `metabolomix/tools/reanalysis_tools.py`（`dataset_mass_accuracy` `dataset_novel_candidates`）
- Test: `tests/test_mass_accuracy.py`

**Interfaces:**
- Produces:
  - `mass_accuracy.sme_ppm(ds, rank=1) -> list[tuple[feature_id, theoretical_mz, ppm]]`（`feature_candidates` の `exp_mass_to_charge` と `theoretical_mass_to_charge`。ppm = (実測 − 理論)/理論 × 1e6）。
  - `summarize(values) -> dict`（`n median q1 q3 within_5ppm within_10ppm`。四分位は `statistics.quantiles(n=4)`（exclusive）。n<2 なら n だけ）。
  - ツール `dataset_mass_accuracy(band_edges: list[float] | None = None, rank: int = 1, dataset_handle=None)`：`{overall, bands: [{min_mz, max_mz, ...}], rank}`（`band_edges=[200]` → `<200` と `≥200`）。許容幅の決定はしない（再解析の規則）。
  - ツール `dataset_novel_candidates(spec: dict, dataset_handle=None, output_path=None)`。spec:

```json
{"statistics": {"result_id": "res_..."} | {"analysis_spec": {<analysis-spec.v1>}, "analysis_id": "main"},
 "evidence": ["sme"],
 "ppm": {"default": 10.0, "bands": [{"max_mz": 200.0, "ppm": 10.0}]},
 "q_max": 0.05,
 "exclude": {"feature_ids": [], "inchikey_blocks": []}}
```
  行: `feature_id name inchikey adduct mz rt ppm group_ref group_test log2fc_test_over_ref p q direction`、q の昇順・feature_id の順。`rule_hash` を返す。

- [ ] **Step 1:** テスト: 手で作った SME の ppm から帯ごとの数字 / 除外と ppm 帯 / result_id と analysis_spec の両方の入口 / 同じ spec で同じ rule_hash。
- [ ] **Step 2–4:** 失敗 → 実装 → 通る。
- [ ] **Step 5:** コミット `feat(analysis): m/z 帯ごとの質量精度と新しい知見の候補の一覧`

### Task 11: 細かい不具合（.d のキャッシュ・メソッドファイルの ASCII）

**Files:**
- Modify: `metabolomix/pipeline/inputs.py`（`_iter_dir_raw_files`）, `metabolomix/console/execution.py`（入力の比較）, `metabolomix/console/method_file.py`, メソッドファイルを計画で読む入口（`console_tools.console_plan` と `pipeline/inputs.py` の計画）
- Test: `tests/test_raw_inventory_bruker_cache.py`, `tests/test_method_file_ascii.py`

**Interfaces:**
- Produces:
  - `inputs.BRUKER_CACHE_NAMES = frozenset({"analysis.sqlite", "analysis.sqlite-journal"})`：`.d` ディレクトリの直下にあるこれらを列挙から除く。記録済みの指紋（変更前の run）との比較では、記録側からも同じ名前を除いてから比べる（古い run が「入力が変わった」にならない）。
  - `method_file.check_ascii(path: Path) -> None`：ASCII 以外のバイトがあれば `DomainError("METHOD_FILE_NOT_ASCII", "...{行}行目 {列}文字目に ASCII 以外の文字 {repr} があります（MS-DIAL Console は ASCII しか読めません）", {"path", "line", "column", "char"})`。最初の 1 件と件数を返す。計画の段階（メソッドファイルを読む最初の関数）で呼ぶ。

- [ ] **Step 1:** テスト: `.d` の中の `analysis.sqlite` があっても無くても指紋が同じ / 記録に sqlite を含む古い指紋と新しい列挙が一致と判定される / 日本語のコメントのメソッドファイルで `console_plan` が `METHOD_FILE_NOT_ASCII`（行・列つき）を返し `'ascii' codec` の例外にならない。
- [ ] **Step 2–4:** 失敗 → 実装 → 通る。
- [ ] **Step 5:** コミット `fix(console): .d の baf2sql キャッシュを入力の指紋から外し、メソッドファイルの ASCII を計画の前に検査する`

### Task 12: 文書と腐敗防止テスト

**Files:**
- Create: `docs/output_format/reanalysis.md`（spec・結果の JSON の欄、KNN/LOESS の参照（R `impute` の版と `impute.knn` の引数、statsmodels の版と `lowess` の引数）、証拠の束、合否、同族列）
- Modify: `docs/output_format/{core,mztab,library,dcl}.md`、`metabolomix/tools/resources.py`（新トピック `reanalysis` の配信）、`docs/workflow/`（新ツール 13 本と、ハンドルの章を `dataset.md` へ。`index.md` の内訳）、`USAGE.md`（13 本と既存の opt-in 引数。件数 90）、`CLAUDE.md` 冒頭の「ツール 77」→ 90 と構成の表に新モジュール 1 行ずつ、`metabolomix/core/mcp_core.py` の `MCP_INSTRUCTIONS` に「並行の再解析はハンドル」を 2〜3 行
- Test: 既存の `tests/test_workflow_docs.py` `tests/test_readme_links.py` `tests/test_server_registration.py`

- [ ] **Step 1:** 文書を直す。
- [ ] **Step 2:** `pytest tests -q` が全部通る。
- [ ] **Step 3:** main ツリーの `docs/HISTRY.md` に追記、`docs/task.md` の該当項目を更新。vault の流れ図（`ms-data-parser-flow.md`）にハンドルの入口と再解析の枝を足し、出典行の日付を直す。
- [ ] **Step 4:** コミット `docs: ハンドルと再解析のツールの文書・出力の定義・流れ`

### 最後に

- [ ] 全体のレビュー（superpowers:requesting-code-review）。
- [ ] main へのマージ（`--no-ff`、`Merge feat/reanalysis-gaps: <要約>`）と push は**ユーザーに確かめてから**。
- [ ] 親のセッションへ、手順書 v7 の実装に要る情報（ツール名・spec の JSON・証拠の束の JSON・ハンドルの使い方）を短く報告する。共有ノートの folder にも同じ内容を置く。
