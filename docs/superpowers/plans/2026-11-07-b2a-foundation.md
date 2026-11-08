# B-2a 土台（描画のロック・読み込みの非同期化・性能・400 検体・1〜2 試料・設定の一本化・据え置き細目）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 公開面（ツールの名前と引数）を変えずに、B-2b の切り替えが乗る土台を直す: ワーカースレッドの描画を直列にし、データセットの読み込みをイベントループから外し、解析フォルダの読み込みを 1 回の読み・ストリーミングにして 400 検体と 1〜2 試料を読めるようにし、設定を `lipidmix.local.toml` に一本化し、B-1 で据え置いた細目を片付ける。

**Architecture:** 描画は `metabolomix/plots/render.py` の 1 つのロックの下でだけ図を作って閉じる。読み込みは `handle_dispatch.ALWAYS_THREADED` でワーカーへ回し、`opening.open_dataset` に入口ごとの読み込み中ロックを置く。`.arf` は `arf/reader.py` のストリーミングの反復子で 1 度だけ読み、`msdial_adapter.load_msdial_folder` は行から直接 numpy に詰める（セルの dict を作らない）。ファイルの sha256 は stdlib だけの新しい leaf `metabolomix/core/file_digest.py` が (パス, 大きさ, mtime) で覚える。設定は `metabolomix/core/user_config.py` の型つきの表 1 つにし、資産と置き場所の環境変数を読まない。

**Tech Stack:** Python 3.14・numpy・msgpack・lz4・matplotlib（Agg）・anyio・pytest。MCP は FastMCP。

**Spec:**
- `docs/superpowers/specs/2026-11-06-unified-dataset-and-tools-design.md`（B。この計画は §7.4 と、`docs/task.md` §2 の B-2 の項目のうち公開面を変えないもの。§8〜§10 と `data_dir` の探索の削除は B-2b `docs/superpowers/plans/2026-11-07-b2b-tool-surface.md`）
- `docs/superpowers/specs/2026-11-06-overhaul-overview-design.md`（全体）

## 決めたこと（設計が任せた点・利用者の判断）

1. **分け方（利用者の判断 2026-11-07）**: B-2 を B-2a（この計画。公開面を変えない）と B-2b（公開面の切り替え）に分け、同じブランチ `claude/overhaul-b2-cdc6e9` で a → b の順に実装する。
2. **置き場所の環境変数（利用者の判断 2026-11-07）**: `LIPIDMIX_REPORTS_DIR` → `[paths] reports_dir`、`LIPIDMIX_LIBRARY_CACHE_DIR` → `[paths] library_cache`、`LIPIDMIX_MAX_DATASET_HANDLES` → `[server] max_dataset_handles`、`LIPIDMIX_CAVEAT_MODE` → `[server] caveat_mode`。`LIPIDMIX_DATA_DIR` は B-2b で探索ごと消す（この計画では触らない）。
3. **資産の環境変数は読まない**: spec §7.4 の 9 つを `user_config` から消す。環境変数が設定ファイルに勝つ規則も消える。テストは設定ファイルを tmp に書いて `LIPIDMIX_CONFIG` で指す。
4. **`[msdial] lbm` が設定ファイルに残っていたら `CONFIG_INVALID`**（「`[library] lbm` へ移した」と書く）。黙って無視すると Console の LBM の解決が変わるため。
5. **ライブラリのキャッシュの既定**は「設定ファイルのフォルダ」の `data/.library-cache`（`user_config.config_file_path().parent / "data" / ".library-cache"`。設定ファイルが無いときも同じ式で、既定の場所なら `<repo>/data/.library-cache`）。今は `mcp_core.DATA_DIR` に付いて動く。worktree は `.mcp.json` の `LIPIDMIX_CONFIG` で main ツリーの設定ファイルを指すので、キャッシュも main ツリーのものを共有する（spec §7.4「worktree は main ツリーの蓄積状態を向く」を、置き場所の環境変数なしで保つ）。
6. **`LIPIDMIX_PIPELINE_INDEX_DIR` は残す**（テストだけが使う上書き。spec §7.4 の一覧に無く、別プロセスの worker の索引をテストごとに分けるため）。文書には書かない。
7. **400 検体**: 読む前に断るのでなく、ストリーミングにして読めるようにする。それでも足りないときだけ、伸長後の大きさと物理メモリから見積もって `ARF_TOO_LARGE`（理由・規模・対処つき）で止める。`MemoryError` も同じ封筒にする（今は `attach_to_dataset` が `arf_unreadable` に飲み込む）。
8. **1〜2 試料**: `_convert_arf_peak_group` の `len(group) >= 3` をやめ、ピークの行の形（`_is_alignment_peak_row`）で群を見分ける。
9. **`dataset_pca` の群の指定（task.md の HOLD）**は C の `dataset_pca`（C §10）へ移す。**分子種ごとの図**は B-2b の `dataset_plot(kind="species" / "species_pca")`。どちらもこの計画では作らない（HISTRY に書く）。
10. **既存の統計関数（`analysis/differential.py` など）は動かさない・消さない**（C のセッションが `metabolomix/analysis/` に足す。衝突を避ける）。

## Global Constraints

- コード・コメント・docstring・文書・コミットメッセージは日本語。
- テストは worktree の中から main ツリーの `.venv` で回す: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。
- コミットは pre-commit が全テストを回す（数分）。`git commit` は最初からバックグラウンドで実行し、出力を scratchpad のファイルへ落として `FAILED` を引く。実行中にファイルを編集しない。`git stash` を使わない。push・`overhaul/integration` へのマージ・main の操作をしない。
- 公開面（ツールの名前・引数・戻り値の欄）を変えない。変えるのは B-2b。例外は設定の案内の文言（環境変数の名前を書いていた箇所）と、新しいエラーコード `ARF_TOO_LARGE`。
- MessagePack の Key 番号は `docs/schema/*.md` から取る。推測で直さない。
- `metabolomix/core/{mcp_core,user_config,atomic_io,process_control,file_digest}.py` は stdlib だけ（`file_digest` は新設の leaf）。`metabolomix/pipeline/{engine,worker,service,recovery,store}.py` から session を import しない。
- 研究室の参照データ（ライブラリ）のパス・ファイル名を追跡対象のテスト・文書に書かない。
- fixture はテスト自身が tmp に作る（`tests/curation_fixtures.py` の `arf_row` `write_arf` `arf2_spot_raw` `write_arf2`）。
- ツールの戻り値は `json_payload()`、全ツール `structured_output=False`、float は丸める。
- 解析フォルダのアダプタの `DatasetState` の中身（値・検出状態・試料・同定・指紋の対象）は、性能の直しの前後で同じ（既存のテスト `tests/test_dataset_agreement.py` `tests/test_folder_tool_parity.py` などがそのまま緑）。
- 新しい環境変数を足さない。

## Review Focus

1. **2 つのデータセットを別々のハンドルで同時に描く**（volcano と mirror を 2 スレッドで並べる）→ どちらの PNG も壊れず、`plt.get_fignums()` が空に戻る。（Task 1 `test_concurrent_renders_are_serialized`）
2. **同じフォルダを 2 つの呼び出しが同時に開く** → 読むのは 1 回、2 つ目は基のキャッシュを使う。別のフォルダどうしは並ぶ。（Task 2 `test_concurrent_open_same_input_reads_once`）
3. **読み込み中に元のファイルが書き換わった** → 基をキャッシュしない（次の読み込みは読み直す）。（Task 3 `test_file_changed_during_load_is_not_cached`）
4. **1 試料・2 試料のアラインメント** → スポットを読み、`.arf2` との数の照合が通る。数が本当に違うときの文言は「別のバッチ」と決めつけない。（Task 6 `test_one_and_two_sample_alignments_read`、`test_spot_count_mismatch_wording`）
5. **設定ファイルに古い `[msdial] lbm` だけが残っている** → `CONFIG_INVALID` で `[library] lbm` へ移すよう案内する。黙って LBM なしで進まない。（Task 8 `test_old_msdial_lbm_key_is_rejected`）

---

### Task 1: 描画のロック

**Files:**
- Modify: `metabolomix/plots/render.py`、`metabolomix/plots/mirror.py`（呼び出しごとの `matplotlib.use("Agg")` を消す）、`metabolomix/plots/result_output.py`、`metabolomix/tools/reports.py`（`_save_eic` `_write_png_and_svg`）、PNG を作る全ツールの呼び出し箇所（`grep -rn "figure_to_png\|savefig" metabolomix`）
- Test: `tests/test_plot_lock.py`

**Interfaces:**
- Produces: `render.PLOT_LOCK: threading.RLock`、`render.render_png(build, /, *args, dpi=DEFAULT_IMAGE_DPI, **kwargs) -> bytes`（ロックの下で `fig = build(*args, **kwargs)` → PNG → `plt.close(fig)`）、`render.render_files(build, paths: dict[str, Path], /, *args, **kwargs) -> None`（同じく、拡張子ごとに `savefig`。`save_figure` 用）。`figure_to_png` はモジュールの中だけで使う（外から呼ばない）。

- [ ] **Step 1: 失敗するテストを書く**: `test_concurrent_renders_are_serialized`（8 スレッドで `render_png(render_volcano_plot, payload)` と `render_png(render_mirror, ...)` を並べ、どの戻り値も `b"\x89PNG"` で始まり、終了後 `plt.get_fignums() == []`。描画関数に `threading.get_ident()` と入退出を記録する包みを当て、入退出が重ならないことを主張）、`test_no_direct_figure_to_png_outside_render`（`metabolomix/` の `.py` を走査し、`render.py` 以外に `figure_to_png(` と `.savefig(` が無い。`arf/reader.py` の CLI 用 `plot_*` は除外リストに名前で載せる）。
- [ ] **Step 2: 走らせて落ちることを確かめる**: `... -m pytest tests/test_plot_lock.py -v` → FAIL（`render_png` が無い）。
- [ ] **Step 3: 実装**: 全呼び出し箇所を `render_png` / `render_files` へ置き換える。`mirror.render_mirror` は内部で `fig.canvas.draw()` を呼ぶので、`render_png(render_mirror, ...)` の形で丸ごとロックの下に入れる。`RLock` にするのは、描画関数の中から別の描画の部品を呼ぶ箇所があってもデッドロックしないため。
- [ ] **Step 4: 走らせて通ることを確かめる**: `tests/test_plot_lock.py` と既存の描画のテスト（`tests/test_volcano_plot.py` `tests/test_plots_mirror.py` `tests/test_eic_multi_plot.py` `tests/test_save_figure.py`）が PASS。
- [ ] **Step 5: コミット**: `fix(plots): ワーカースレッドの描画を 1 つのロックで直列にする`

### Task 2: データセットの読み込みをイベントループから外す

**Files:**
- Modify: `metabolomix/core/handle_dispatch.py`（`ALWAYS_THREADED` に `dataset_load`）、`metabolomix/dataset/opening.py`（入口ごとの読み込み中ロック）、`metabolomix/core/session_state.py`（セッションのライブラリを移す箇所を `_current_lock` の下に）、`metabolomix/library/tools.py`（`library_load` がセッションのライブラリを書く箇所を同じロックの下に）
- Test: `tests/test_load_off_event_loop.py`

**Interfaces:**
- Produces: `opening._inflight_lock(key: str) -> threading.Lock`（基のキャッシュのキーごと。辞書は `threading.Lock` で守る）。`open_dataset` はキーのロックを取ってから `cached_base` を見直し、無ければ読む。

- [ ] **Step 1: 失敗するテストを書く**:
  - `test_dataset_load_does_not_block_loop`: `opening._open` を 0.5 秒眠る偽物に差し替え、`asyncio` で `mcp.call_tool("dataset_load", ...)` と、0.05 秒ごとに数を足す協程を並べる。読み込みの間に協程が 5 回以上進む。
  - `test_concurrent_open_same_input_reads_once`: 同じ `directory` を 2 スレッドで `open_dataset`、`_open` の呼び出しは 1 回、2 つ目の `reused_base` が真。別の 2 つのフォルダでは `_open` が並ぶ（両方の入時刻が両方の出時刻より前）。
  - `test_library_move_and_library_load_do_not_race`: 読み込みの setter がライブラリを移す間に `library_load` を呼んでも、セッションのライブラリがちょうど 1 つのハンドルに属する（どちらの順でも）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装。** 旧 ARF 経路の `load_dataset`（`tools/dataset.py`）はワーカーへ回さない（`session.arf` と `mcp_core.DATA_DIR` をロックなしで書き、イベントループの ARF のツールが同時に読むとちぎれるため。B-2b で新しい `load_dataset` に置き換わり、そちらを `ALWAYS_THREADED` にする）。
- [ ] **Step 4: 通ることを確かめる**（`tests/test_session_handle_unification.py` `tests/test_dataset_handles.py` を含む）。
- [ ] **Step 5: コミット**: `fix(dataset): 読み込みをワーカーへ回し、同じ入口の同時の読み込みを 1 回にする`

### Task 3: ファイルの sha256 の記憶と、指紋を読む前に取る

**Files:**
- Create: `metabolomix/core/file_digest.py`
- Modify: `metabolomix/dataset/library_memory.py`（`_measure`）、`metabolomix/mztab/dataset_state.py`（`_sha256`）、`metabolomix/handoff/schema.py`（`sha256_file`）、`metabolomix/library/store.py`（`source_sha256`。永続の `.library-digests.json` は残し、メモリの記憶を先に引く）、`metabolomix/dataset/opening.py`（`open_dataset`）
- Test: `tests/test_file_digest.py`、`tests/test_dataset_opening.py`（追記）

**Interfaces:**
- Produces: `file_digest.sha256_of(path: str | Path) -> str`（キー `(resolve した文字列, st_size, st_mtime_ns)`、`threading.Lock` で守る、上限 256 件の LRU）、`file_digest.clear() -> None`（テスト用）。stdlib だけ。

- [ ] **Step 1: 失敗するテストを書く**: `test_same_file_hashed_once`（`builtins.open` を数える包みで 2 回目は開かない）、`test_changed_size_or_mtime_rehashes`、`test_file_changed_during_load_is_not_cached`（`_open` の途中でフォルダの `.arf2` に 1 バイト足す偽物にし、`registry.cached_base(key)` が None。次の `open_dataset` の `reused_base` が偽）、`test_core_file_digest_is_stdlib_only`（`tests/test_user_config.py` の import 検査と同じ方法）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装。** `open_dataset` は `_open` の前に `_dependencies` の stat を取り、後の stat と違えば `put_base` しない（そのまま返す。`notes` に「読み込みの間に入力が変わった」）。B-1 の据え置き「指紋を読み込み後に取る」を閉じる。`_dependencies` は読む前には分からないファイル（`.mddata` の後出しなど）があるので、読む前の stat は「読んだ後の `_dependencies` の集合のうち、読む前に stat できたもの」で比べる。
- [ ] **Step 4: 通ることを確かめる**（`tests/test_library_memory.py` `tests/test_library_store.py` `tests/test_dataset_opening.py`）。
- [ ] **Step 5: コミット**: `perf(dataset): ファイルの sha256 を (パス・大きさ・時刻) で覚え、読み込み中の変更をキャッシュしない`

### Task 4: `.arf` を 1 回だけ読み、セルの dict を作らない

**Files:**
- Modify: `metabolomix/arf/reader.py`（ストリーミングの反復子）、`metabolomix/mztab/evidence.py`（`load_arf_evidence` `matching_arfs` `normalize_arf_spots` の利用者）、`metabolomix/dataset/opening.py`（`_matching_arf` `_folder_via_job`）、`metabolomix/dataset/msdial_adapter.py`（`_read_groups` `_cells` `_check_alignment` `load_msdial_folder`）
- Test: `tests/test_arf_streaming.py`、`tests/test_dataset_opening.py`（追記）

**Interfaces:**
- Produces:
  - `arf_reader.iter_groups(file_like) -> Iterator[dict]`: `deserialize` と同じ要素（`{"MasterAlignmentID", "MassCenter", "PeakProperties": [slim row, ...]}` など今の `deserialize` の 1 要素の形）を 1 つずつ返す。msgpack の `Unpacker` で大きな配列の要素を 1 つずつ取り出し、全体のリストを作らない。`deserialize(f)` は `list(iter_groups(f))` にする（既存の呼び出しは変えない）。
  - `arf_reader.uncompressed_size(file_like) -> int`: LZ4 の各ブロックの先頭の整数（伸長後の大きさ）の和。伸長しない。
  - `evidence.ArfScan`（dataclass）: `path`、`spot_mz: list[float | None]`（`spot_representative_mz` の定義）、`detected: np.ndarray`、`values: dict[str, np.ndarray]`（measure ごと）、`per_sample_mz`、`file_ids: list[int]`、`file_names: dict[int, str | None]`、`msms_assigned: np.ndarray`。`.arf` を 1 回なめて作る。
  - `evidence.scan_arf(path, *, measures: tuple[str, ...]) -> ArfScan`（`arf_reader.peak_field` で行から直接取り、`alignment_feature_row` を呼ばない）。
  - `load_msdial_folder(selection, *, measure="peak_height", scan: ArfScan | None = None)`。`scan` が渡れば `.arf` を読まない。
  - `_matching_arf(...)` は合格した候補の `ArfScan` を返す（今の戻り値の組に 1 つ足す）。2 本目が合格したら `ALIGNMENT_AMBIGUOUS` なので、持つのは最初の合格だけ。
- [ ] **Step 1: 失敗するテストを書く**: `test_iter_groups_matches_deserialize`（`write_arf` の fixture で `list(iter_groups) == deserialize`）、`test_iter_groups_is_lazy`（`inspect.isgenerator`、かつ 1 つ目を取った時点で最後の群を変換していない: 変換関数を数える包み）、`test_folder_load_builds_no_cell_dicts`（`arf_reader.alignment_feature_row` を例外を投げる偽物にして `load_msdial_folder` が通る）、`test_job_folder_reads_arf_once`（`arf_reader.iter_groups` を数える包みで job の解析フォルダの経路が 1 回）、`test_scan_equals_old_cells`（既存の fixture 群で `scan_arf` の値・検出・m/z が今の `_cells` の結果と一致）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装。** `_read_groups` の `io.BytesIO(handle.read())` をやめてファイルを直に渡す。`_check_alignment` は `ArfScan.spot_mz` と `.arf2` の m/z を比べる（今と同じ判定 `compare_feature_mz`）。`cells_by_spot` を消す。`attach_to_dataset`（mzTab-M のアダプタ）も `scan_arf` を使う。
- [ ] **Step 4: 通ることを確かめる**（`tests/test_dataset_agreement.py` `tests/test_folder_tool_parity.py` `tests/test_mztab_evidence.py` ほか全体）。
- [ ] **Step 5: コミット**: `perf(dataset): .arf を 1 回のストリーミングで読み、セルの dict を作らない`

### Task 5: 400 検体（大きさの見積もりと `ARF_TOO_LARGE`）と基のキャッシュの上限

**Files:**
- Modify: `metabolomix/core/mcp_errors.py`（エラーコード）、`metabolomix/mztab/evidence.py`（`scan_arf` の前の見積もり、`attach_to_dataset` の `except` を `MemoryError` で分ける）、`metabolomix/dataset/msdial_adapter.py`、`metabolomix/core/handles.py`（`MAX_BASES` と大きさ）、`metabolomix/arf/reader.py`（`_ARF_BYTES_PER_ROW` を今の行の大きさに直すか、使わなくなれば消す）
- Test: `tests/test_arf_scale.py`、`tests/test_dataset_handles.py`（追記）

**Interfaces:**
- Produces:
  - `evidence.estimate_scan_bytes(path) -> int`（`uncompressed_size` ＋ 行あたりの固定費 × 見積もりの行数。行数は伸長後の大きさ ÷ 660 バイト。定数は根拠をコメントに書く: HISTRY 2026-10-15 の 966 MB / 1.47M 行）
  - `core.handles.available_memory_bytes() -> int | None`（Windows は `ctypes` の `GlobalMemoryStatusEx` の `ullAvailPhys`、POSIX は `os.sysconf("SC_AVPHYS_PAGES") * SC_PAGE_SIZE`、取れなければ None）。`handles` に置くのは基のキャッシュも使うため。stdlib だけ。
  - エラー `ARF_TOO_LARGE`（`DomainError`）: `details = {"path_name", "estimated_bytes", "available_bytes", "remedy"}`。`remedy` は「試料を分けてアラインメントし直す・メモリの多い機械で読む・mzTab-M のアダプタ（`source="mztab"`）なら値は読めるが検出状態は付かない」。見積もりが空きメモリの 80% を超えたら読む前に投げる。読んでいる途中の `MemoryError` も同じコードにする。
  - `handles.MAX_BASE_BYTES = 2 * 1024**3`、`handles.estimate_state_bytes(ds) -> int`（`feature_matrix` `detected_mask` `pp_matrix` と `analysis_matrices` の配列の `nbytes` ＋ `per_sample_mz` の要素数 × 100 ＋ `feature_metadata` の件数 × 1500）。`put_base` は件数 4 と合計バイトの両方で LRU を追い出し、1 つで上限を超える基はキャッシュしない。
- [ ] **Step 1: 失敗するテストを書く**: `test_too_large_refused_before_read`（`available_memory_bytes` を小さな値に差し替え、`load_msdial_folder` が `ARF_TOO_LARGE` を投げ、`iter_groups` を呼んでいない）、`test_memory_error_becomes_too_large`（`iter_groups` の途中で `MemoryError`）、`test_mztab_attach_reports_too_large_not_unreadable`（mzTab-M のアダプタの注記が `arf_too_large`）、`test_400_samples_stream_within_bound`（合成の 300 スポット × 400 試料の `.arf` を `tracemalloc` の下で読み、ピークが行列の大きさ `300*400*9` バイトの 40 倍未満、かつ値が `deserialize` 経由の旧い組み立てと一致）、`test_base_cache_size_cap`（`MAX_BASE_BYTES` を小さく差し替え、大きい基が追い出される・上限を超える基は入らない）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装。** 旧 ARF 経路の `ArfMemoryError` と `_format_arf_memory_error`（`arf/tools.py`）はこの計画では残す（B-2b で消える）。
- [ ] **Step 4: 通ることを確かめる。**
- [ ] **Step 5: コミット**: `feat(dataset): 400 検体の .arf をストリーミングで読み、読めない大きさは ARF_TOO_LARGE で返す`

### Task 6: 1〜2 試料のアラインメントと、スポット数の食い違いの文言

**Files:**
- Modify: `metabolomix/arf/reader.py`（`_convert_arf_peak_group` `extract_arf_data` `_is_arf_peak_group`。使われない `_is_arf_feature_container` を消す）、`metabolomix/dataset/msdial_adapter.py`（`_check_alignment` の文言）、`metabolomix/mztab/evidence.py`（`build_evidence` の `feature_count_mismatch` の注記）、`tests/test_assay_evidence.py`（制限を書いた箇所）
- Test: `tests/test_arf_small_alignments.py`

- [ ] **Step 1: 失敗するテストを書く**: `test_one_and_two_sample_alignments_read`（`write_arf` で 1 試料・2 試料の fixture を作り、`deserialize` がスポットの数だけ群を返し、`load_msdial_folder` が `n_samples == 1 / 2` で開く）、`test_header_lists_are_not_groups`（ChromXs の形の入れ子リスト `[[1.0,0,0],[2.0,0,0],...]` を先頭に置いた fixture で群に数えない）、`test_spot_count_mismatch_wording`（数が本当に違うとき、`ALIGNMENT_MISMATCH` の文言が「.arf のスポット数（N）と .arf2 のスポット数（M）が違う。同じアラインメントの出力でないか、.arf を読み切れていない。」で、「別のバッチ」と決めつけない）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装。** 群の条件は「空でないリストで、要素がすべて `_is_alignment_peak_row`（`len(row) > 25` かつ Key18 が数）」にする（`docs/schema/AlignmentChromPeakFeature.md` で Key18 を確かめる）。
- [ ] **Step 4: 通ることを確かめる**（`tests/test_parser_decode.py` `tests/test_arf_hot_paths.py` を含む全体）。
- [ ] **Step 5: コミット**: `fix(arf): 1〜2 試料のアラインメントを読む`

### Task 7: mzTab-M のアダプタの `.arf2` の見つけ方

**Files:**
- Modify: `metabolomix/curation/apply.py`（`arf2_for_mztab`）
- Test: `tests/test_curation_apply.py`（追記）

**Interfaces:**
- Produces: `arf2_for_mztab(mztab_path, *, feature_mz: list[float | None] | None = None) -> Path | None`。今の「ファイル名にバッチの語幹」の一致が 1 本ならそれ。無ければ同じフォルダ（と job の `spot_catalog`）の `.arf2` のうち、スポット数が SMF の数と同じで `compare_feature_mz` が通るものがちょうど 1 本ならそれ。それ以外は None（今と同じ）。呼び出し側（`tools/dataset_analysis_tools.py` の 3 か所と `tools/handle_tools.alignment_files_for`）は `feature_mz` を `ds.feature_metadata` の m/z から渡す。
- [ ] **Step 1: 失敗するテストを書く**: `test_arf2_found_without_batch_stem`（mzTab の名前に語幹が無い fixture で、フォルダの `.arf2` が 1 本で m/z が合えば見つかり、`curation.state` が `None` でない）、`test_arf2_not_guessed_when_mz_differs`。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る**（`tests/test_folder_tool_parity.py` の mzTab の名前の回避は外してよい）。
- [ ] **Step 5: コミット**: `fix(curation): mzTab-M の名前にバッチの語幹が無くても .arf2 を m/z の照合で見つける`

### Task 8: 設定の型と、資産の環境変数の削除・LBM の設定の一本化

**Files:**
- Modify: `metabolomix/core/user_config.py`、`metabolomix/console/runner.py`（`get_exe_path`）、`metabolomix/tools/console_tools.py`（`_lbm_setting` `_msp_setting` `_msdial_exe_setup_help` と案内の文言）、`metabolomix/pipeline/inputs.py`（`_resolve_lbm_pinned`）、`metabolomix/core/path_resolvers.py`（`LIBRARY_ENV_VARS` を設定のキーの表へ）、`metabolomix/pathway/bridge.py`、`metabolomix/console/method_file.py` `metabolomix/library/{tools,store}.py` `metabolomix/tools/handle_tools.py` の文言、`lipidmix.example.toml`、`tests/conftest.py`、`tests/pipeline_fixtures.py`、環境変数を使うテスト（`grep -rln "MSDIAL_EXE\|MSDIAL_LBM\|MSDIAL_MSP_\|MSDIAL_LIBRARY_LBM\|MSLIPIDMAPPER_" tests`）
- Create: `tests/config_fixtures.py`
- Test: `tests/test_user_config.py`（書き直し）、`tests/test_no_asset_env_vars.py`

**Interfaces:**
- Produces:
  - `user_config.SETTINGS: dict[str, str]`（キー → 種類 `"path"` / `"text"` / `"int"`）。キー: `msdial.exe` `library.msp_positive` `library.msp_negative` `library.lbm` `mslipidmapper.rscript` `mslipidmapper.browser`（path）、`mslipidmapper.mslipidmapper_api_url` `mslipidmapper.api_token`（text）。Task 9 が `paths.*` `server.*` を足す。
  - `MOVED_KEYS = {"msdial.lbm": "library.lbm"}`: 設定ファイルにあれば `ConfigInvalidError`（文言「[msdial] lbm は [library] lbm へ移しました」）。
  - `Setting(key, value: str | int, config_file)`（`source` と `env_var` を消す）。`setting_label` は「<file> の [sec] name」だけ、`missing_hint` / `describe_missing` は設定ファイルの書き方だけを案内する（`details` の `env_var` を消し `config_key` を足す）。
  - `get_setting(key) -> Setting | None`: 設定ファイルだけを読む。種類 `int` は TOML の整数（bool は拒む）、`path` は今どおり設定ファイルのフォルダから解く、`text` は前後の空白を落とす。
  - `tests/config_fixtures.py`: `use_config(monkeypatch, tmp_path, values: dict[str, object]) -> Path`（ドットのキーを節に分けて TOML を書き、パスはリテラル文字列 `'...'`、`LIPIDMIX_CONFIG` を向け、`user_config` のキャッシュを消す）。
- [ ] **Step 1: 失敗するテストを書く**: `test_env_var_is_ignored`（`MSDIAL_EXE` を置いても `get_setting("msdial.exe")` は設定ファイルの値か None）、`test_old_msdial_lbm_key_is_rejected`、`test_int_setting_rejects_bool_and_string`、`test_console_lbm_reads_library_lbm`（`[library] lbm` だけを書いた設定で `console_plan` の LBM の解決が 5 番目の段で設定を使う。解決の順は変わらない）、`test_no_asset_env_vars`（`metabolomix/` の `.py` に 9 つの名前の文字列が `os.environ` / `getenv` / 案内の文言として現れない。エラーコード `MSDIAL_EXE_NOT_FOUND` などは名前が違うので当たらない正規表現 `\b(MSDIAL_EXE|MSDIAL_LBM|MSDIAL_MSP_POS|MSDIAL_MSP_NEG|MSDIAL_LIBRARY_LBM|MSLIPIDMAPPER_[A-Z_]+)\b(?!_NOT)` で、`MSDIAL_EXE_` に続くものを除く）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装**し、環境変数を `setenv` していたテストを `use_config` へ移す（主張は変えない。`test_console_method_file.py` の `Setting(key="msdial.lbm", ...)` は `library.lbm` へ）。`conftest._isolate_external_asset_settings` は `delenv` の輪をやめ、`LIPIDMIX_CONFIG` を存在しないパスへ向けることだけを残す。`_msdial_exe_setup_help` の PowerShell の `SetEnvironmentVariable` の手順を、`lipidmix.local.toml` の `[msdial] exe` の書き方に替える（再起動は要らない、と書く）。
- [ ] **Step 4: 通ることを確かめる**（全体）。
- [ ] **Step 5: コミット**: `refactor(config): 外部資産の場所は設定ファイルだけで指し、LBM の設定を [library] lbm にまとめる`

### Task 9: 置き場所の環境変数を設定ファイルへ

**Files:**
- Modify: `metabolomix/core/user_config.py`（キーを足す）、`metabolomix/core/mcp_core.py`（`_report_dir_candidates`）、`metabolomix/library/store.py`（`cache_dir`、`LIBRARY_CACHE_ENV` を消す）、`metabolomix/core/handles.py`（`max_handles`、`MAX_HANDLES_ENV` を消す。`DATASET_HANDLE_LIMIT` の文言）、`metabolomix/core/session_state.py`（`maybe_prepend_caveat`）、`lipidmix.example.toml`、これらの環境変数を使うテスト（`REPORTS_DIR` 14 ファイル、`LIBRARY_CACHE` 約 22 ファイル、`MAX_DATASET_HANDLES` 2、`CAVEAT_MODE` 2）、`tests/curation_fixtures.py`
- Test: `tests/test_user_config.py`（追記）

**Interfaces:**
- Produces: `SETTINGS` に `paths.reports_dir` `paths.library_cache`（path）、`server.max_dataset_handles`（int、1 未満は 1）、`server.caveat_mode`（text、`digest` / `off`。他の値は `CONFIG_INVALID`）。`store.cache_dir()` の既定は `user_config.config_file_path().parent / "data" / ".library-cache"`（決めたこと 5）。`tests/conftest.py` は `LIPIDMIX_CONFIG` をテストごとの tmp の中の存在しないパスへ向ける（既定のキャッシュがリポジトリの中にできないように）。`_report_dir_candidates()` は `[reports_dir, DATA_DIR/reports]`（`reports_dir` が無ければ後者だけ。意味は今のまま。B-2b で解析フォルダへ替わる）。
- [ ] **Step 1: 失敗するテストを書く**: `test_place_settings_from_config`（4 つのキーがそれぞれの利用者に効く）、`test_place_env_vars_are_ignored`（4 つの環境変数を置いても効かない）、`test_library_cache_does_not_follow_data_dir`（`mcp_core.DATA_DIR` を差し替えてもキャッシュの場所は変わらない）、`test_library_cache_default_follows_config_folder`（`LIPIDMIX_CONFIG` を別のフォルダの設定ファイルへ向けると、既定のキャッシュがそのフォルダの `data/.library-cache`）、`test_caveat_mode_invalid_value_rejected`。`tests/test_no_asset_env_vars.py` の正規表現に `LIPIDMIX_(REPORTS_DIR|LIBRARY_CACHE_DIR|MAX_DATASET_HANDLES|CAVEAT_MODE)` を足す。
- [ ] **Step 2〜4: 落ちる → 実装（テストを `use_config` へ移す）→ 通る。**
- [ ] **Step 5: 文書**: `CLAUDE.md` の環境変数の段落（残る環境変数は `LIPIDMIX_CONFIG` `LIPIDMIX_DATA_DIR`（B-2b まで）`LIPIDMIX_TRANSPORT` `LIPIDMIX_HOST` `LIPIDMIX_PORT` `LIPIDMIX_PLOT_OUTPUT`、`[library] lbm` の一本化、環境変数が優先の記述を消す）、`README.md` `USAGE.md` `DEPLOY.md` `docs/cli.md` `docs/output_format/{library,reanalysis}.md` `docs/workflow/{library,pipeline,plots}.md` `.gitignore` のコメント、`lipidmix.example.toml` の注を直す。`docs/superpowers/` の過去の設計書・計画書は履歴なので直さない。
- [ ] **Step 6: コミット**: `refactor(config): 置き場所の環境変数を lipidmix.local.toml の [paths] [server] へ移す`

### Task 10: B-1 の据え置き細目（ファイルの選択・同定の欄・試料・注記）

**Files:**
- Modify: `metabolomix/core/path_resolvers.py`（`_pick_latest`）、`metabolomix/msdial/alignment_spot.py`（`alignment_sibling`）、`metabolomix/dataset/opening.py`（`_folder_sidecars`、`pipeline_arf2`、`_BROKEN_MZTAB_HINT`）、`metabolomix/arf2/spot_fields.py`、`metabolomix/arf2/match_results.py`、`metabolomix/dataset/msdial_adapter.py`、`metabolomix/dataset/msdial_identity.py`、`metabolomix/dataset/library_memory.py`（`is_lbm_annotator`）、`metabolomix/analysis/dataset_analysis.py`（`_DATE_RE` と `import re`）、`metabolomix/analysis/claim_evidence.py` か証拠の束の組み立て（`msms_not_evaluable_reason` の位置）、`tests/test_msdial_classes.py`
- Test: 下の各テストを該当するテストファイルに足す

各項目と、閉じる条件（テスト名）:

| 項目 | 直し方 | テスト |
|---|---|---|
| 時刻の無い `.arf` が `_pick_latest` で勝つ | 時刻を持つ候補があれば、時刻の無い候補より先に選ぶ（PeakProperties の優先はその後で当てる）。時刻の無い候補どうしは mtime | `test_timestamped_batch_beats_untimestamped` |
| `.arf2` の探索と `alignment_sibling` の大文字小文字 | `alignment_sibling` と `_folder_sidecars` もフォルダを走査し、名前を `casefold()` で比べる | `test_sidecar_found_case_insensitively` |
| `.arf2` 入力で PeakProperties / DriftSpots が並ぶ | テストだけ足す（今の規則で PeakProperties） | `test_arf2_input_prefers_peakproperties` |
| 空のスポット配列の `int(raw_spot[0])` | `spot_fields` と `match_results` で空・短い配列を `DomainError("ARF2_SPOT_MALFORMED")`（スポットの位置つき） | `test_empty_spot_array_is_reported` |
| `spot_fields` の値の主張 | テストだけ足す（Key ごとの値を `docs/schema/AlignmentSpotProperty.md` に照らして主張） | `test_spot_fields_values` |
| `internal_standard_id` の -1 | **直さない**（B-3 の定量で意味を決める。task.md に残す） | — |
| `test_msdial_classes.py` の mddata の包みの重複 | 共通の fixture 関数にまとめる | 既存のテストが緑 |
| `per_sample_mz` の空 dict と非正の m/z | m/z が 0 以下・NaN は持たない。検出セルが 1 つも無い特徴は空 dict を持たない（キーごと省く） | `test_per_sample_mz_skips_nonpositive` |
| FileID の改名の注記が空白だけの違いで出る | 比べる前に両方を `strip()` | `test_rename_note_ignores_whitespace` |
| SML の (d) が空白だけの名前を通す | 名前を `strip()` してから空・`Unknown`・`null` を判定 | `test_whitespace_name_is_not_sml` |
| `identity_columns` の "unknown" の大文字小文字 | `casefold()` で比べる | `test_unknown_name_case_insensitive` |
| `is_lbm_annotator` が False を返さない | 末尾 `\.msp\d*$`（大文字小文字を問わない）なら False、`.lbm\d*` なら True、ほかは None | `test_is_lbm_annotator_false_for_msp` |
| pipeline の退避が最新の `input/` のバッチしか見ない | `pipeline_arf2` を新しい順の候補列 `pipeline_arf2_candidates(pipeline_path) -> list[Path]` にし、`_matching_arf` は合うまで順に試す（二重に数えない R28 は保つ） | `test_pipeline_falls_back_to_older_input_batch` |
| 壊れた mzTab-M の案内 | `_BROKEN_MZTAB_HINT` に、読めなかった理由（例外の要約）と `source="mztab"` でも読めないことを足す | `test_broken_mztab_hint_names_reason` |
| `_DATE_RE` の残り | `dataset_analysis.py` の `_DATE_RE` と `import re` を消す（`sample_manifest.py` が正準） | 既存のテストが緑 |
| 証拠の束で `msms_not_evaluable_reason` が先頭のキー | 束の dict の末尾へ置く | `test_not_evaluable_reason_is_last_key` |
| 閉じたハンドルの競合で `required_tools: dataset_load` | B-2b の改名（`load_dataset`）で直す。この計画では触らない | — |
| `session.arf` をワーカーから読むときのちぎれ | B-2b でスロットごと消える | — |

- [ ] **Step 1: 失敗するテストを書く**（上の表のテスト）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装。**
- [ ] **Step 4: 通ることを確かめる**（全体）。
- [ ] **Step 5: コミット**: `fix(dataset): B-1 の据え置き細目（選択・大文字小文字・同定の欄・注記）`

### Task 11: 記録

**Files:**
- Modify（main ツリー側。追跡外）: `C:\Users\yuu18\Metabolomix_with_LLM\docs\HISTRY.md`（末尾に「2026-11-07 B-2a 完了」節）、`C:\Users\yuu18\Metabolomix_with_LLM\docs\task.md`（末尾への追記だけ。並列中なので既存行の書き換え・削除をしない）

- [ ] **Step 1: HISTRY**: 決めたこと（この計画の冒頭 1〜10）、実装中の裁定、性能の前後（解析フォルダの読み込みの所要。fixture の数字でなく、手元に実データがあれば `MTBLS9365/pipeline/neg/a3777451` を `open_dataset(pipeline_path, source="folder")` で開いた秒を前後で 1 度ずつ。研究データに書き込まないことを mtime で確かめる）、**流れが変わる点**（vault の流れ図は親セッションが直すので列挙だけ: 資産と置き場所の環境変数が効かなくなる、`[msdial] lbm` の拒否、`ARF_TOO_LARGE`、1〜2 試料が読める、mzTab-M の `.arf2` の見つけ方）。
- [ ] **Step 2: task.md の末尾に追記**: 「B-2a で閉じた項目」（task.md §2 の該当行の見出しを引用して、閉じたことを書く。行そのものは並列の作業が終わってから整理する）と、残り（`internal_standard_id` の -1 は B-3、`dataset_pca` の群は C）。
- [ ] **Step 3: 全テスト**: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q` の結果を HISTRY に書く。
