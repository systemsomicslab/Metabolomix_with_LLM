# library_mode="msp_only" Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `console_plan`・`console_method_template`・pipeline v1 に `library_mode="msp_only"` を足し、LBM を使わず研究室 MSP（極性で pos / neg）だけで MS-DIAL Console を回せるようにする。

**Architecture:** 純粋関数（`method_file.resolve_msp` / `msp_only_overrides`）が MSP の解決と実効メソッドの上書きを決め、Console 系ツールと pipeline v1 の `inspect_inputs` がそれを呼ぶ。既定の `library_mode="auto"` は今の LBM 解決の経路をそのまま通る。研究室 MSP の置き場所（ディレクトリ）は戻り値・エラー文に載せず、ファイル名だけを返す。

**Tech Stack:** Python 3.14（`C:/Python314/python.exe`）、pytest、FastMCP（既存）。

**Spec:** `docs/superpowers/specs/2026-10-07-library-mode-msp-only-design.md`

## Global Constraints

- 作業場所は worktree `C:\Users\yuu18\Metabolomix_with_LLM\.worktrees\feat-library-mode-msp-only`（ブランチ `feat/library-mode-msp-only`）。テストはそのルートから `C:/Python314/python.exe -m pytest <path> -q`。
- **コミットのたびに pre-commit が全テスト（約 3,000 件・約 4 分）を走らせる。** `git commit` は最初からバックグラウンドで実行し、出力をファイルに落として `FAILED` を後から引く。実行中にファイルを編集しない（フックは作業ツリーを検証する）。`--no-verify` は使わない。
- `git stash` を使わない（worktree 間で共有される）。
- 研究室ライブラリの実パス・実ファイル名をテスト・追跡対象の文書に書かない（R5）。テストは tmp に偽の `.msp` を作り、設定は tmp の toml を `LIPIDMIX_CONFIG` で指す（`tests/conftest.py` の autouse fixture が `MSDIAL_*` を消し、`LIPIDMIX_CONFIG` を存在しないパスへ向けている）。
- **研究室 MSP の置き場所（ディレクトリ）を戻り値・エラー文・実行レポートに載せない。** ファイル名・設定キー・環境変数名・設定ファイルのパスだけを載せる（`core/path_resolvers.py` `LibraryPathError` と同じ規約）。絶対パスは実効メソッドとデータ側の入力記録（`inputs` / snapshot）にだけ置く。
- `library_mode="auto"`（既定）の振る舞い・戻り値の既存キー・要求の内容 hash（`tests/test_metabolomics_stages.py` `_V1_DEFAULT_FINGERPRINT`）を変えない。
- 文書は日本語、識別子は英語。コミットメッセージの末尾に `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`。
- 記録は main ツリーの `C:\Users\yuu18\Metabolomix_with_LLM\docs\HISTRY.md`（追記）と `docs\task.md`（Task 9 の項目を消す）。マージは `--no-ff`、`Merge feat/library-mode-msp-only: <日本語の要約>`。push はユーザーに都度確認する。vault の流れ図は直さない（引数だけの変更）。

## Review Focus

1. **NAS 上の MSP が一時的に読めない**（sha256 計算中の `OSError`）— pipeline_plan は封筒（`MSP_NOT_FOUND`）で止まるべきで、例外で落ちてはいけない。Task 4 でテストする。
2. **メソッドファイルに `Msp file path` や `Text DB file path` が重複行・小文字の綴りで書かれている** — 全行が差し替わるべき（Console は後勝ちで読む）。Task 1 でテストする。
3. **極性の値が大文字始まり**（`"Negative"`）で計画に入る — 設定キーの引き当てで `KeyError` にならず、小文字にして引くべき。Task 4 でテストする。
4. **`omics="metabolomics"` で `msp_only`** — リピドミクスと同じく LBM が空・MSP が入るべき。Task 2 でテストする。
5. **環境変数 `MSDIAL_MSP_*` が設定されている** — 設定ファイルより優先され、`source` が `"env"` になるべき（契約上は置かない値なので、出どころが見えることが大事）。Task 2 でテストする。

## File Structure

| ファイル | 責務 | Task |
|---|---|---|
| `metabolomix/core/user_config.py` | 極性 → MSP 設定キーの正準 `MSP_SETTING_KEYS` | 1 |
| `metabolomix/core/path_resolvers.py` | `LIBRARY_SETTING_KEYS` を `MSP_SETTING_KEYS` の参照に変える | 1 |
| `metabolomix/console/method_file.py` | `MspResolution` / `resolve_msp` / `IDENTIFICATION_CLEAR_KEYS` / `msp_only_overrides` | 1 |
| `metabolomix/tools/console_tools.py` | `console_plan` / `console_method_template` の `library_mode` / `msp_file` | 2 |
| `tests/lab_msp_fixtures.py`（新規） | tmp に偽の研究室 MSP と設定ファイルを作る helper | 2 |
| `metabolomix/pipeline/request.py` / `request_v2.py` | 要求キー・検証・内容 hash の既定値除外・v2 の拒否 | 3 |
| `metabolomix/pipeline/inputs.py` / `service.py` / `report.py` / `tools/pipeline_tools.py` | 入力計画・検査・指紋・receipt・レポート | 4 |
| `USAGE.md` / `docs/workflow/pipeline.md` | 利用者向け文書 | 2・4 |

---

### Task 1: MSP の解決と実効メソッドの上書き（純粋関数）

**Files:**
- Modify: `metabolomix/core/user_config.py`（`SETTINGS` の直後）
- Modify: `metabolomix/core/path_resolvers.py:253`（`LIBRARY_SETTING_KEYS`）
- Modify: `metabolomix/console/method_file.py`（import、`LbmResolution` の直後、`resolve_lbm` の直後）
- Test: `tests/test_console_method_file.py`（末尾に追加）

**Interfaces:**
- Produces:
  - `user_config.MSP_SETTING_KEYS: dict[str, str]` = `{"positive": "library.msp_positive", "negative": "library.msp_negative"}`
  - `method_file.MspResolution`（frozen dataclass: `path: str | None`, `source: str`, `error_code: str | None = None`, `message: str | None = None`）。`source` は `argument` / `env` / `config_file` / `not_configured`
  - `method_file.resolve_msp(override: str | None, polarity: str, msp_setting: Setting | None) -> MspResolution`
  - `method_file.IDENTIFICATION_CLEAR_KEYS: tuple[str, ...]` = `(TEXT_DB_KEY, *SETTINGS_PATH_KEYS)`
  - `method_file.msp_only_overrides(method_keys: dict[str, str], msp_path: str) -> tuple[dict[str, str], list[str]]`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_console_method_file.py` の末尾に追加:

```python
# ---------- library_mode="msp_only": MSP の解決と実効メソッドの上書き ----------

import os

from metabolomix.console import method_file as _mf
from metabolomix.core import path_resolvers as _path_resolvers
from metabolomix.core import user_config as _user_config


def _msp_setting(path, *, source="config_file", key="library.msp_negative"):
    return _user_config.Setting(
        key=key, value=str(path), source=source,
        env_var=_user_config.SETTINGS[key],
        config_file=None if source == "env" else "C:/cfg/lipidmix.local.toml")


def test_msp_setting_keys_are_the_single_source_for_library_lookup():
    assert _user_config.MSP_SETTING_KEYS == {
        "positive": "library.msp_positive", "negative": "library.msp_negative"}
    assert _path_resolvers.LIBRARY_SETTING_KEYS is _user_config.MSP_SETTING_KEYS


def test_resolve_msp_prefers_the_argument(tmp_path):
    arg = tmp_path / "given.msp"
    arg.write_text("NAME: a\n", encoding="ascii")
    other = tmp_path / "configured.msp"
    other.write_text("NAME: b\n", encoding="ascii")
    res = _mf.resolve_msp(str(arg), "negative", _msp_setting(other))
    assert res == _mf.MspResolution(path=os.path.abspath(arg), source="argument")


def test_resolve_msp_uses_the_setting_and_reports_where_it_came_from(tmp_path):
    lib = tmp_path / "lab_neg.msp"
    lib.write_text("NAME: a\n", encoding="ascii")
    from_file = _mf.resolve_msp(None, "negative", _msp_setting(lib))
    assert from_file.path == os.path.abspath(lib) and from_file.source == "config_file"
    from_env = _mf.resolve_msp(None, "negative", _msp_setting(lib, source="env"))
    assert from_env.source == "env"


def test_resolve_msp_not_configured_names_the_setting_key_and_env_var():
    res = _mf.resolve_msp(None, "positive", None)
    assert res.error_code == "MSP_NOT_CONFIGURED"
    assert res.path is None and res.source == "not_configured"
    assert "msp_positive" in res.message and "MSDIAL_MSP_POS" in res.message


def test_resolve_msp_not_found_shows_only_the_file_name(tmp_path):
    hidden_dir = tmp_path / "secret_lab_share"
    missing = hidden_dir / "lab_neg.msp"
    res = _mf.resolve_msp(None, "negative", _msp_setting(missing))
    assert res.error_code == "MSP_NOT_FOUND"
    assert "lab_neg.msp" in res.message
    assert "secret_lab_share" not in res.message
    arg = _mf.resolve_msp(str(missing), "negative", None)
    assert arg.error_code == "MSP_NOT_FOUND" and arg.source == "argument"
    assert "secret_lab_share" not in arg.message


def test_msp_only_overrides_blank_lbm_set_msp_and_clear_other_identification_keys():
    method_keys = {
        "ion mode": "Negative",
        "lbm file path": "C:/libs/lipids.lbm2",
        "msp file path": "C:/public/other.msp",
        "text db file path": "db.txt",                    # 小文字の綴りでも拾う
        "msp annotator settings file path": "msp.tsv",
        "isotope text db file path": "iso.txt",           # 同定用ではないので残す
        "compounds library file path for rt correction": "rt.txt",
    }
    overrides, removed = _mf.msp_only_overrides(method_keys, "C:/lab/lab_neg.msp")
    assert overrides[_mf.LBM_KEY] == ""
    assert overrides[_mf.MSP_KEY] == "C:/lab/lab_neg.msp"
    assert overrides[_mf.TEXT_DB_KEY] == ""
    assert overrides["MSP annotator settings file path"] == ""
    assert "Isotope text DB file path" not in overrides
    assert _mf.RT_REFERENCE_KEY not in overrides
    assert removed == [_mf.LBM_KEY, _mf.TEXT_DB_KEY, "MSP annotator settings file path"]


def test_msp_only_overrides_does_not_add_keys_the_method_never_declared():
    overrides, removed = _mf.msp_only_overrides({"ion mode": "Positive", "lbm file path": ""},
                                                "C:/lab/lab_pos.msp")
    assert overrides == {_mf.LBM_KEY: "", _mf.MSP_KEY: "C:/lab/lab_pos.msp"}
    assert removed == []


def test_blank_override_replaces_every_duplicate_line(tmp_path):
    src = tmp_path / "params.txt"
    src.write_text(
        "Ion mode: Negative\n"
        "Lbm file path: C:/libs/a.lbm2\n"
        "MSP FILE PATH: C:/public/one.msp\n"
        "Lbm file path: C:/libs/b.lbm2\n"
        "msp file path: C:/public/two.msp\n",
        encoding="ascii")
    overrides, _ = _mf.msp_only_overrides(_mf.read_method_keys(src), "C:/lab/lab_neg.msp")
    dest = _mf.write_effective_method_file(src, tmp_path / "out" / "effective.txt", overrides)
    lines = dest.read_text(encoding="ascii").splitlines()
    assert lines.count("Lbm file path: ") == 2
    assert lines.count("Msp file path: C:/lab/lab_neg.msp") == 2
    assert not any("one.msp" in line or "two.msp" in line or ".lbm2" in line for line in lines)
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_console_method_file.py -q -k "msp"`
Expected: FAIL（`AttributeError: module 'metabolomix.core.user_config' has no attribute 'MSP_SETTING_KEYS'` などで 8 件落ちる）

- [ ] **Step 3: 実装する**

`metabolomix/core/user_config.py` — `SETTINGS` の定義の直後に追加:

```python
#: 極性 → 研究室参照ライブラリ（.msp）の設定キー。正準はここ
#: （`core/path_resolvers.LIBRARY_SETTING_KEYS` と `console/method_file.resolve_msp` が参照する。
#: pipeline の worker も使うので、mcp_core を引き込まない leaf のここに置く）。
MSP_SETTING_KEYS: dict[str, str] = {
    "positive": "library.msp_positive",
    "negative": "library.msp_negative",
}
```

`metabolomix/core/path_resolvers.py:253` を置き換える:

```python
LIBRARY_SETTING_KEYS = user_config.MSP_SETTING_KEYS
```

`metabolomix/console/method_file.py` — import 行を置き換える:

```python
from metabolomix.core.user_config import MSP_SETTING_KEYS, Setting, missing_hint, setting_label
```

`LbmResolution` の直後に追加:

```python
@dataclass(frozen=True)
class MspResolution:
    """library_mode="msp_only" の MSP の解決結果。`error_code` が非 None なら呼び出し側は停止する。

    `message` に研究室ライブラリの置き場所（ディレクトリ）を載せない——戻り値は LLM の
    文脈に入る（`core/path_resolvers.LibraryPathError` と同じ規約）。ファイル名だけを言う。
    """

    path: str | None
    source: str  # argument | env | config_file | not_configured
    error_code: str | None = None
    message: str | None = None
```

`SETTINGS_PATH_KEYS` の定義の直後に追加:

```python
#: library_mode="msp_only" で実効メソッドから消す同定用の宣言。研究室 MSP 以外の
#: ライブラリ（Text DB、複数 MSP / Text の注釈器設定表）が混ざらないようにする。
#: 同定に使わないキー（アイソトープ追跡・ターゲット検出・RT 補正）は含めない。
IDENTIFICATION_CLEAR_KEYS: tuple[str, ...] = (TEXT_DB_KEY, *SETTINGS_PATH_KEYS)
```

`resolve_lbm` の直後に追加:

```python
_MSP_WHAT = "研究室の参照ライブラリ（.msp）"


def resolve_msp(override: str | None, polarity: str, msp_setting: Setting | None) -> MspResolution:
    """library_mode="msp_only" の MSP を決める。順に: 明示引数 → 極性の設定。

    極性の設定は呼び出し側が `user_config.get_setting(MSP_SETTING_KEYS[polarity])` で引いて
    渡す（環境変数 MSDIAL_MSP_POS / MSDIAL_MSP_NEG → lipidmix.local.toml）。メソッドファイルの
    `Msp file path` の宣言は使わない——研究室 MSP 以外を宣言した GUI パラメータを流用しても
    契約から外れないようにする（spec 2026-10-07 §3.2）。返すパスは絶対パス。
    """
    if override:
        candidate = Path(override).expanduser()
        if candidate.is_file():
            return MspResolution(path=_absolute(candidate), source="argument")
        return MspResolution(
            path=None, source="argument", error_code="MSP_NOT_FOUND",
            message=f"msp_file が指すファイル（{candidate.name}）がありません。")
    key = MSP_SETTING_KEYS[polarity]
    if msp_setting is None:
        return MspResolution(
            path=None, source="not_configured", error_code="MSP_NOT_CONFIGURED",
            message=(f"library_mode='msp_only' には {polarity} の{_MSP_WHAT}が要ります。"
                     + missing_hint(key, _MSP_WHAT)))
    candidate = Path(msp_setting.value).expanduser()
    if candidate.is_file():
        return MspResolution(path=_absolute(candidate), source=msp_setting.source)
    label = setting_label(msp_setting)
    return MspResolution(
        path=None, source=msp_setting.source, error_code="MSP_NOT_FOUND",
        message=f"{label} が指すファイル（{candidate.name}）がありません。{label} を確認してください。")


def msp_only_overrides(method_keys: dict[str, str], msp_path: str) -> tuple[dict[str, str], list[str]]:
    """library_mode="msp_only" の実効メソッドの上書きと、消した宣言のキーを返す。

    `Lbm file path` は空、`Msp file path` は `msp_path`。原本が空でなく宣言している
    `IDENTIFICATION_CLEAR_KEYS` は空にする。キーは正準の綴りで返す
    （`write_effective_method_file` は大文字小文字を問わず全行を差し替える）。
    `method_keys` は `read_method_keys` の戻り値（キーは小文字）。
    """
    overrides: dict[str, str] = {LBM_KEY: "", MSP_KEY: msp_path}
    removed: list[str] = []
    if (method_keys.get(LBM_KEY.lower()) or "").strip():
        removed.append(LBM_KEY)
    for key in IDENTIFICATION_CLEAR_KEYS:
        if (method_keys.get(key.lower()) or "").strip():
            overrides[key] = ""
            removed.append(key)
    return overrides, removed
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_console_method_file.py tests/test_user_config.py tests/test_library_store.py -q`
Expected: PASS（全件。`test_library_store.py` が無ければ `tests/test_library*.py` を対象にする）

- [ ] **Step 5: コミットする（バックグラウンド）**

```bash
cd /c/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-library-mode-msp-only
git add metabolomix/core/user_config.py metabolomix/core/path_resolvers.py metabolomix/console/method_file.py tests/test_console_method_file.py
git commit -m "feat(console): library_mode=msp_only 用に MSP の解決と実効メソッドの上書きを足す

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$TMP/commit_task1.log" 2>&1; echo "exit=$?"; grep -E "passed|failed|FAILED" "$TMP/commit_task1.log" | tail -5
```

Expected: `exit=0`、全テスト通過。

---

### Task 2: console_plan / console_method_template に library_mode と msp_file

**Files:**
- Create: `tests/lab_msp_fixtures.py`
- Modify: `metabolomix/tools/console_tools.py`（`console_plan` 23〜322 行、`console_method_template` 662〜783 行、末尾の helper 群）
- Modify: `USAGE.md:136,140`（`console_plan` / `console_method_template` の行）
- Test: `tests/test_console_plan_method.py`（末尾に追加）

**Interfaces:**
- Consumes: `method_file.resolve_msp`、`method_file.msp_only_overrides`、`user_config.MSP_SETTING_KEYS`（Task 1）
- Produces:
  - `console_plan(..., library_mode: str = "auto", msp_file: str | None = None)`、`console_method_template(..., library_mode: str = "auto", msp_file: str | None = None)`
  - 戻り値に `library_mode`、`lbm`、`msp: {"file", "source"}`、`removed_declarations`
  - `tests/lab_msp_fixtures.write_lab_msp_config(tmp_path, monkeypatch, *, positive=True, negative=True, directory_name="lab_library") -> dict[str, Path]`（Task 4 も使う）

- [ ] **Step 1: 共通 fixture を作る**

`tests/lab_msp_fixtures.py`:

```python
"""library_mode="msp_only" のテスト用: tmp に偽の研究室 MSP と設定ファイルを作る。

研究室ライブラリの実パス・実ファイル名は使わない（外部流出禁止の資産）。設定ファイルは
tmp に書いて `LIPIDMIX_CONFIG` で指す（tests/conftest.py の隔離 fixture の流儀）。
"""
from __future__ import annotations

from pathlib import Path


def write_lab_msp_config(tmp_path: Path, monkeypatch, *, positive: bool = True,
                         negative: bool = True, directory_name: str = "lab_library") -> dict[str, Path]:
    lib_dir = tmp_path / directory_name
    lib_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    lines = ["[library]"]
    for polarity, enabled, key in (("positive", positive, "msp_positive"),
                                   ("negative", negative, "msp_negative")):
        if not enabled:
            continue
        path = lib_dir / f"lab_{polarity[:3]}.msp"
        path.write_text("NAME: fake\nPRECURSORMZ: 100\nNum Peaks: 0\n", encoding="ascii")
        lines.append(f"{key} = '{path}'")
        paths[polarity] = path
    config = tmp_path / "lipidmix.local.toml"
    config.write_text("\n".join(lines) + "\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(config))
    return paths
```

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_console_plan_method.py` の末尾に追加:

```python
# ---------- library_mode="msp_only"（LBM 不使用・研究室 MSP のみ） ----------

import os

from tests.lab_msp_fixtures import write_lab_msp_config


def _effective_keys(path) -> dict[str, list[str]]:
    keys: dict[str, list[str]] = {}
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            keys.setdefault(key.strip().lower(), []).append(value.strip())
    return keys


def _msp_only_plan(tmp_path, monkeypatch, method_text, **kwargs):
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")  # auto なら拾われる LBM
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text(method_text, encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    kwargs.setdefault("polarity", "negative")
    kwargs.setdefault("library_mode", "msp_only")
    return _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                    measure="peak_height", **kwargs))


def test_console_plan_msp_only_blanks_lbm_and_writes_the_lab_msp(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(
        tmp_path, monkeypatch,
        "Ion mode: Negative\nLbm file path: \nMsp file path: C:/public/other.msp\n"
        "Text DB file path: other.txt\n")
    assert parsed["status"] == "planned"
    assert parsed["library_mode"] == "msp_only"
    assert parsed["lbm"] == {"path": None, "source": "disabled"}
    assert parsed["msp"] == {"file": "lab_neg.msp", "source": "config_file"}
    assert parsed["removed_declarations"] == ["Text DB file path"]
    keys = _effective_keys(parsed["method_file"])
    assert keys["lbm file path"] == [""]
    assert keys["msp file path"] == [os.path.abspath(paths["negative"])]
    assert keys["text db file path"] == [""]
    assert "lab_library" not in _json.dumps(parsed, ensure_ascii=False)  # 置き場所を返さない


def test_console_plan_msp_only_picks_the_positive_library_for_positive(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Positive\n", polarity="positive")
    assert parsed["msp"]["file"] == "lab_pos.msp"
    assert _effective_keys(parsed["method_file"])["msp file path"] == [os.path.abspath(paths["positive"])]


def test_console_plan_msp_only_for_metabolomics(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\nLbm file path: C:/x/y.lbm2\n",
                            omics="metabolomics")
    assert parsed["status"] == "planned"
    assert parsed["removed_declarations"] == ["Lbm file path"]
    assert _effective_keys(parsed["method_file"])["lbm file path"] == [""]


def test_console_plan_msp_only_env_var_wins_and_is_reported(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    env_lib = tmp_path / "env_dir" / "env_neg.msp"
    env_lib.parent.mkdir()
    env_lib.write_text("NAME: e\n", encoding="ascii")
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(env_lib))
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["msp"] == {"file": "env_neg.msp", "source": "env"}


def test_console_plan_msp_only_argument_wins(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    given = tmp_path / "given.msp"
    given.write_text("NAME: g\n", encoding="ascii")
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n", msp_file=str(given))
    assert parsed["msp"] == {"file": "given.msp", "source": "argument"}


def test_console_plan_msp_only_stops_when_not_configured(tmp_path, monkeypatch):
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["error"]["code"] == "MSP_NOT_CONFIGURED"
    assert not (tmp_path / "runs").exists()  # ジョブを作らない


def test_console_plan_msp_only_stops_when_the_file_is_missing(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    paths["negative"].unlink()
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["error"]["code"] == "MSP_NOT_FOUND"
    assert "secret_share" not in _json.dumps(parsed, ensure_ascii=False)
    assert not (tmp_path / "runs").exists()


def test_console_plan_msp_only_rejects_lbm_file(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n", lbm_file="C:/x/y.lbm2")
    assert parsed["error"]["code"] == "LIBRARY_MODE_CONFLICT"


def test_console_plan_auto_rejects_msp_file(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch)
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n",
                            library_mode="auto", msp_file="C:/x/y.msp")
    assert parsed["error"]["code"] == "LIBRARY_MODE_CONFLICT"


def test_console_plan_rejects_an_unknown_library_mode(tmp_path, monkeypatch):
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n", library_mode="lbm_only")
    assert parsed["error"]["code"] == "JOB_NOT_PLANNED"


def test_console_plan_msp_only_rejects_a_non_ascii_library_path(tmp_path, monkeypatch):
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="研究室")
    parsed = _msp_only_plan(tmp_path, monkeypatch, "Ion mode: Negative\n")
    assert parsed["error"]["code"] == "METHOD_ENCODING_UNSUPPORTED"
    assert "研究室" not in _json.dumps(parsed, ensure_ascii=False)
    assert not (tmp_path / "runs").exists()


def test_console_plan_auto_reports_msp_not_used(tmp_path, monkeypatch):
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\nLbm file path: \n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["library_mode"] == "auto"
    assert parsed["lbm"]["source"] == "exe_dir"
    assert parsed["msp"] == {"file": None, "source": "not_used"}
    assert parsed["removed_declarations"] == []


def test_console_method_template_msp_only(tmp_path, monkeypatch):
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    exe = _fake_exe_with_lbm(tmp_path, "Msp_lipids.lbm2")
    _plan_ready(tmp_path, monkeypatch, exe)
    src = tmp_path / "neg_param.txt"
    src.write_text("Ion mode: Negative\nLbm file path: C:/x/y.lbm2\nText DB file path: db.txt\n",
                   encoding="ascii")
    out = tmp_path / "pos_param.txt"
    from metabolomix.tools.console_tools import console_method_template
    parsed = _json.loads(console_method_template(
        out_path=str(out), polarity="positive", based_on=str(src), library_mode="msp_only"))
    assert parsed["status"] == "written"
    assert parsed["lbm"] == {"path": None, "source": "disabled"}
    assert parsed["msp"] == {"file": "lab_pos.msp", "source": "config_file"}
    assert parsed["removed_declarations"] == ["Lbm file path", "Text DB file path"]
    keys = _effective_keys(out)
    assert keys["lbm file path"] == [""]
    assert keys["msp file path"] == [os.path.abspath(paths["positive"])]
    assert keys["ion mode"] == ["Positive"]
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_console_plan_method.py -q -k "msp or library_mode"`
Expected: FAIL（`TypeError: console_plan() got an unexpected keyword argument 'library_mode'` など）

- [ ] **Step 4: 実装する**

`metabolomix/tools/console_tools.py` — 末尾の helper 群（`_lbm_setting` の直後）に追加:

```python
_LIBRARY_MODES = ("auto", "msp_only")


def _library_mode_error(library_mode: str, lbm_file: str | None, msp_file: str | None) -> str | None:
    """library_mode の値と、lbm_file / msp_file との組み合わせを検査する。"""
    if library_mode not in _LIBRARY_MODES:
        return console_error("JOB_NOT_PLANNED",
                             f"library_mode は 'auto' または 'msp_only' です: {library_mode!r}")
    if library_mode == "msp_only" and lbm_file:
        return console_error(
            "LIBRARY_MODE_CONFLICT",
            "library_mode='msp_only' は LBM を使わないので、lbm_file と同時には指定できません。",
            {"library_mode": library_mode})
    if library_mode == "auto" and msp_file:
        return console_error(
            "LIBRARY_MODE_CONFLICT",
            "msp_file は library_mode='msp_only' のときだけ使えます。library_mode='auto' では"
            "MSP を解決しないので、渡しても使われません。",
            {"library_mode": library_mode})
    return None


def _msp_setting(polarity: str):
    """極性の研究室 MSP の設定。設定ファイルが読めなければ封筒（str）を返す。"""
    from metabolomix.core import user_config
    try:
        return user_config.get_setting(user_config.MSP_SETTING_KEYS[polarity])
    except user_config.ConfigInvalidError as exc:
        return console_error(exc.code, exc.message, exc.details())


def _resolve_libraries(method_keys: dict, method_path: Path, *, omics: str, exe: str | None,
                       polarity: str, library_mode: str, lbm_file: str | None = None,
                       msp_file: str | None = None):
    """実効メソッドの上書き（相対宣言の絶対化＋ライブラリ）と戻り値の断片を返す。

    失敗したら封筒（str）。研究室 MSP の置き場所は戻り値・エラー文に載せない
    （ファイル名だけ。core/path_resolvers.LibraryPathError と同じ規約）。
    """
    from metabolomix.console import method_file as method_file_mod
    from metabolomix.core import user_config

    overrides = method_file_mod.relative_path_overrides(method_keys, method_path)
    if library_mode == "msp_only":
        msp_setting = _msp_setting(polarity)
        if isinstance(msp_setting, str):
            return msp_setting
        msp = method_file_mod.resolve_msp(msp_file, polarity, msp_setting)
        if msp.error_code:
            details = ({"msp_file": Path(msp_file).name} if msp.source == "argument"
                       else user_config.describe_missing(
                           user_config.MSP_SETTING_KEYS[polarity], msp_setting))
            return console_error(msp.error_code, msp.message or "", details)
        if not msp.path.isascii():
            return console_error(
                "METHOD_ENCODING_UNSUPPORTED",
                f"研究室 MSP（{Path(msp.path).name}）のパスに ASCII 以外の文字があります。"
                "MS-DIAL Console はメソッドファイルを ASCII で読むため、このライブラリは"
                "見つからず黙って同定 0 件になります。ASCII だけのフォルダへ置き直してください。",
                {"keys": [method_file_mod.MSP_KEY]})
        extra, removed = method_file_mod.msp_only_overrides(method_keys, msp.path)
        overrides.update(extra)
        return overrides, {
            "library_mode": "msp_only",
            "lbm": {"path": None, "source": "disabled"},
            "msp": {"file": Path(msp.path).name, "source": msp.source},
            "removed_declarations": removed,
        }

    lbm_setting = _lbm_setting()
    if isinstance(lbm_setting, str):
        return lbm_setting
    lbm = method_file_mod.resolve_lbm(
        method_keys, method_path, omics=omics, exe_path=exe, lbm_setting=lbm_setting,
        override=lbm_file)
    if lbm.error_code:
        return console_error(lbm.error_code, lbm.message or "",
                             {"candidates": list(lbm.candidates)} if lbm.candidates else None)
    if lbm.path:
        overrides[method_file_mod.LBM_KEY] = lbm.path
    return overrides, {
        "library_mode": "auto",
        "lbm": {"path": lbm.path, "source": lbm.source},
        "msp": {"file": None, "source": "not_used"},
        "removed_declarations": [],
    }
```

`console_plan`:
- シグネチャの `keep_extension: str | None = None,` の後に `library_mode: str = "auto",` と `msp_file: str | None = None,` を足す。
- docstring の `keep_extension:` の節の後に追加:

```
    library_mode:
        "auto"（既定）か "msp_only"。"msp_only" は LBM を使わず、研究室の参照ライブラリ
        （.msp）だけで同定する: 実効メソッドの `Lbm file path` を空にし、`Msp file path` に
        `msp_file` か極性の設定（環境変数 MSDIAL_MSP_POS / MSDIAL_MSP_NEG →
        lipidmix.local.toml の [library] msp_positive / msp_negative）を書き、メソッドが宣言する
        Text DB と MSP / Text の注釈器設定表を空にする（消したキーは removed_declarations）。
        メソッドの `Msp file path` の宣言は使わない。MSP が決まらなければ MSP_NOT_CONFIGURED /
        MSP_NOT_FOUND で止まり、ジョブを作らない。lbm_file とは同時に指定できない。
    msp_file:
        library_mode="msp_only" で使う MSP の明示パス（設定より優先）。"auto" では指定できない。
```

- `if omics not in (...)` の直後に追加:

```python
    library_error = _library_mode_error(library_mode, lbm_file, msp_file)
    if library_error:
        return library_error
```

- 244〜259 行（`lbm_setting = _lbm_setting()` から `overrides[method_file_mod.LBM_KEY] = lbm.path` まで）を次に置き換える:

```python
    resolved_libraries = _resolve_libraries(
        method_keys, mf, omics=omics, exe=exe, polarity=polarity, library_mode=library_mode,
        lbm_file=lbm_file, msp_file=msp_file)
    if isinstance(resolved_libraries, str):
        return resolved_libraries
    overrides, library_payload = resolved_libraries
```

（直前のコメント「解決した LBM と、相対で宣言されたライブラリ等のパスは…」は残し、`overrides = method_file_mod.relative_path_overrides(...)` の行は消す——`_resolve_libraries` が作る。）
- 戻り値の `"lbm": {"path": lbm.path, "source": lbm.source},` を `**library_payload,` に置き換える。

`console_method_template`:
- シグネチャの `omics: str = "lipidomics",` の後に `library_mode: str = "auto",` と `msp_file: str | None = None,` を足し、docstring に `console_plan` と同じ 2 項を足す（「解決順は console_plan と同じ」と 1 行で書いてよい）。
- `if polarity not in ...` の検査の直後に追加:

```python
    library_error = _library_mode_error(library_mode, None, msp_file)
    if library_error:
        return library_error
```

- `overrides = {**method_file_mod.relative_path_overrides(src_keys, src), ION..., ADDUCT...}` を
  `overrides = {ION_MODE_KEY: ..., ADDUCT_KEY: ...}` の形に分け、`lbm_setting = _lbm_setting()` から
  `overrides[method_file_mod.LBM_KEY] = lbm.path` までを次に置き換える:

```python
    resolved_libraries = _resolve_libraries(
        src_keys, src, omics=omics, exe=exe, polarity=polarity, library_mode=library_mode,
        msp_file=msp_file)
    if isinstance(resolved_libraries, str):
        return resolved_libraries
    library_overrides, library_payload = resolved_libraries
    overrides = {
        **library_overrides,
        method_file_mod.ION_MODE_KEY: polarity.capitalize(),
        method_file_mod.ADDUCT_KEY: method_file_mod.STANDARD_ADDUCTS[polarity],
    }
```

- 戻り値の `"lbm": {"path": lbm.path, "source": lbm.source},` を `**library_payload,` に置き換える。

`USAGE.md:136`（`console_plan` の行）: 引数一覧に `library_mode`, `msp_file` を足し、行末の直前に次の文を足す:
「**`library_mode="msp_only"` は LBM を使わず研究室の参照ライブラリ(.msp)だけで同定する** — 実効メソッドの `Lbm file path` を空にし、`Msp file path` に `msp_file` か極性の設定(`MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`)を書き、メソッドが宣言する Text DB と MSP / Text の注釈器設定表を空にする(戻り値 `removed_declarations`)。メソッドの `Msp file path` の宣言は使わない。決まらなければ `MSP_NOT_CONFIGURED` / `MSP_NOT_FOUND` でジョブを作らずに止まる。`lbm_file` との同時指定と、`auto` での `msp_file` は `LIBRARY_MODE_CONFLICT`。戻り値の `msp` はファイル名と出どころだけで、ライブラリの置き場所は返さない。」
`USAGE.md:140`（`console_method_template` の行）: 引数一覧に `library_mode`, `msp_file` を足し、「`library_mode="msp_only"` は `console_plan` と同じく LBM を空にし研究室 MSP を書く。」を足す。

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_console_plan_method.py tests/test_console_method_template.py tests/test_readme_links.py tests/test_server_registration.py -q`
Expected: PASS（全件）

- [ ] **Step 6: コミットする（バックグラウンド）**

```bash
cd /c/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-library-mode-msp-only
git add metabolomix/tools/console_tools.py tests/lab_msp_fixtures.py tests/test_console_plan_method.py USAGE.md
git commit -m "feat(console): console_plan と console_method_template に library_mode=msp_only を足す

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$TMP/commit_task2.log" 2>&1; echo "exit=$?"; grep -E "passed|failed|FAILED" "$TMP/commit_task2.log" | tail -5
```

Expected: `exit=0`、全テスト通過。

---

### Task 3: pipeline-request.v1 の library_mode / msp_file

**Files:**
- Modify: `metabolomix/pipeline/request.py`（`_TOP_LEVEL_KEYS` 78〜82 行、`_NULL_REJECTED_TOP_LEVEL_KEYS` 138〜141 行、`_validate_fields` 323 行〜、`resolve_request` の `data = {...}` 610〜624 行、`request_fingerprint` 696〜700 行）
- Modify: `metabolomix/pipeline/request_v2.py:636-640`
- Modify: `tests/test_pipeline_request.py`（`test_defaults_are_fully_populated` に 2 行足す・末尾に追加）

**Interfaces:**
- Produces: 解決済みの要求に `library_mode: "auto" | "msp_only"`（既定 `"auto"`）と `msp_file: str | None`（既定 None）。Task 4 の `inspect_inputs` が `request["library_mode"]` / `request["msp_file"]` を読む。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_pipeline_request.py` の `test_defaults_are_fully_populated` の `assert req["lbm_file"] is None` の直後に追加:

```python
    assert req["library_mode"] == "auto"
    assert req["msp_file"] is None
```

末尾に追加:

```python
# ---------- library_mode / msp_file（spec 2026-10-07） ----------

def _root(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    return root


def test_library_mode_msp_only_is_accepted_with_an_msp_file(tmp_path):
    req = resolve_request(_root(tmp_path), {"library_mode": "msp_only", "msp_file": "lab.msp"})
    assert req["library_mode"] == "msp_only"
    assert req["msp_file"] == "lab.msp"
    assert req["value_sources"]["library_mode"] == "explicit"


@pytest.mark.parametrize("explicit", [
    {"library_mode": "lbm_only"},
    {"library_mode": None},
    {"msp_file": None},
    {"library_mode": "msp_only", "msp_file": ""},
    {"library_mode": "msp_only", "lbm_file": "x.lbm2"},
    {"msp_file": "lab.msp"},  # auto では使われないので黙って受けない
])
def test_invalid_library_mode_combinations_are_rejected(tmp_path, explicit):
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve_request(_root(tmp_path), explicit)


def test_library_mode_cannot_be_changed_on_resume(tmp_path):
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(resolve_request(_root(tmp_path)), {"library_mode": "msp_only"})


def test_default_library_mode_does_not_change_the_request_fingerprint(tmp_path):
    from metabolomix.core.atomic_io import canonical_hash
    from metabolomix.pipeline import request as request_mod
    req = resolve_request(_root(tmp_path))
    legacy_keys = request_mod._TOP_LEVEL_KEYS - {"library_mode", "msp_file"}
    assert request_fingerprint(req) == canonical_hash({k: req[k] for k in legacy_keys if k in req})


def test_msp_only_changes_the_request_fingerprint(tmp_path):
    root = _root(tmp_path)
    assert (request_fingerprint(resolve_request(root, {"library_mode": "msp_only"}))
            != request_fingerprint(resolve_request(root)))


def test_v2_rejects_library_mode_with_a_migration_hint(tmp_path):
    from metabolomix.pipeline import request_v2
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID") as info:
        request_v2._check_known_keys({"library_mode": "msp_only"})
    assert "library_mode" in str(info.value)
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pipeline_request.py -q`
Expected: FAIL（`KeyError: 'library_mode'`、未知キーとしての拒否など）

- [ ] **Step 3: 実装する**

`metabolomix/pipeline/request.py`:

`_TOP_LEVEL_KEYS` に `"library_mode", "msp_file"` を足す:

```python
_TOP_LEVEL_KEYS = frozenset({
    "schema", "target", "method_file", "lbm_file", "polarity", "measure",
    "keep_extension", "timeout_s", "save_project", "output_root",
    "sample_manifest", "preprocess", "comparisons",
    # spec 2026-10-07（library_mode="msp_only"）
    "library_mode", "msp_file",
})
```

`_TARGET_VALUES` の直後に追加:

```python
_LIBRARY_MODE_VALUES = frozenset({"auto", "msp_only"})

#: 内容 hash から外す既定値（spec 2026-10-07 §3.1）。既定値のキーが増えただけで既存の
#: 要求の hash が変わると、過去の run への再送が IDEMPOTENCY_CONFLICT になり、
#: 完了済み run の再利用も外れる（tests/test_metabolomics_stages.py が実測値で縛る）。
_FINGERPRINT_OMITTED_DEFAULTS = {"library_mode": "auto", "msp_file": None}
```

`_NULL_REJECTED_TOP_LEVEL_KEYS` に `"library_mode", "msp_file"` を足す。

`_validate_fields` の `for field in ("method_file", "lbm_file", "output_root", "sample_manifest"):` を
`for field in ("method_file", "lbm_file", "msp_file", "output_root", "sample_manifest"):` に変え、その直後に追加:

```python
    library_mode = data.get("library_mode")
    if library_mode not in _LIBRARY_MODE_VALUES:
        _fail(f"library_modeが不正です: {library_mode!r}", library_mode=library_mode)
    if library_mode == "msp_only" and data.get("lbm_file"):
        _fail("library_mode='msp_only'はLBMを使わないので、lbm_fileと同時には指定できません。",
              library_mode=library_mode)
    if library_mode == "auto" and data.get("msp_file"):
        _fail("msp_fileはlibrary_mode='msp_only'のときだけ使えます（autoではMSPを解決しない）。",
              library_mode=library_mode)
```

`resolve_request` の `data = {...}` の `"lbm_file": _pick("lbm_file"),` の直後に追加:

```python
        "library_mode": _pick("library_mode", "auto"),
        "msp_file": _pick("msp_file"),
```

`request_fingerprint` の最後の 2 行を置き換える:

```python
    content = {key: request[key] for key in key_set if key in request}
    if key_set is _TOP_LEVEL_KEYS:
        for key, default in _FINGERPRINT_OMITTED_DEFAULTS.items():
            if key in content and content[key] == default:
                del content[key]
    return canonical_hash(content)
```

（`UPDATABLE` は変えない。`library_mode` / `msp_file` は resume で変えると `NEW_PIPELINE_REQUIRED` になる既存の経路に乗る——Step 2 の `test_library_mode_cannot_be_changed_on_resume` がそれを確かめる。乗らずに別のエラーになったら、`merge_updates` の UPDATABLE 外キーの扱いを読んで同じ経路に乗せる。）

`metabolomix/pipeline/request_v2.py:636-640` を置き換える:

```python
    if {"method_file", "lbm_file", "library_mode", "msp_file"} & unknown:
        _fail(
            "pipeline-request.v2ではmethod_file/lbm_file/library_mode/msp_fileを直接指定"
            "できません（profile_fileから解決してください）。", unknown_keys=sorted(unknown),
        )
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pipeline_request.py tests/test_metabolomics_stages.py tests/test_metabolomics_request.py -q`
Expected: PASS（全件。`test_v1_request_fingerprint_is_unchanged` も変更なしで通る）

- [ ] **Step 5: コミットする（バックグラウンド）**

```bash
cd /c/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-library-mode-msp-only
git add metabolomix/pipeline/request.py metabolomix/pipeline/request_v2.py tests/test_pipeline_request.py
git commit -m "feat(pipeline): pipeline-request.v1 に library_mode と msp_file を足す（既定値は内容 hash に入れない）

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$TMP/commit_task3.log" 2>&1; echo "exit=$?"; grep -E "passed|failed|FAILED" "$TMP/commit_task3.log" | tail -5
```

Expected: `exit=0`、全テスト通過。

---

### Task 4: pipeline v1 の入力計画・検査・receipt・レポート

**Files:**
- Modify: `metabolomix/pipeline/inputs.py`（`_resolve_lbm_pinned` の直後に `_resolve_msp_pinned`、`inspect_inputs` 457〜522 行、`verify_inputs` 782〜785 行）
- Modify: `metabolomix/pipeline/service.py`（`_plan_fingerprint` 115〜127 行、`_resolved_settings` 180〜199 行、`plan_pipeline` の docstring）
- Modify: `metabolomix/pipeline/report.py`（`_section_method` 388〜400 行）
- Modify: `metabolomix/tools/pipeline_tools.py`（`pipeline_plan` docstring の `request:` と `resolved` の説明）
- Modify: `USAGE.md:199`（`pipeline_plan` の行）、`docs/workflow/pipeline.md`（入力計画の節）
- Test: `tests/test_pipeline_inputs.py`、`tests/test_pipeline_service.py`（既存の `test_plan_receipt_shows_the_resolved_settings_the_instructions_promise` の集合を更新・追加）、`tests/test_pipeline_report.py`

**Interfaces:**
- Consumes: `method_file.resolve_msp` / `msp_only_overrides`（Task 1）、`user_config.MSP_SETTING_KEYS`（Task 1）、`request["library_mode"]` / `request["msp_file"]`（Task 3）、`tests.lab_msp_fixtures.write_lab_msp_config`（Task 2）
- Produces: 入力計画（と snapshot）に `library_mode`、`msp: {"path", "file", "sha256", "source"}`、`removed_declarations`。receipt の `resolved` に `library_mode` と `msp: {"file", "sha256", "source"}`。

- [ ] **Step 1: 失敗するテストを書く（入力計画と検査）**

`tests/test_pipeline_inputs.py` の末尾に追加:

```python
# ---------- library_mode="msp_only"（spec 2026-10-07） ----------

from tests.lab_msp_fixtures import write_lab_msp_config


def _msp_only_plan(tmp_path, monkeypatch, **explicit):
    _allow_fake_exe(monkeypatch)
    paths = write_lab_msp_config(tmp_path, monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"], {"library_mode": "msp_only", **explicit})
    return src, paths, inspect_inputs(src["root"], request, exe_path=src["exe"])


def test_msp_only_plan_pins_the_lab_msp_and_blanks_the_lbm(tmp_path, monkeypatch):
    src, paths, plan = _msp_only_plan(tmp_path, monkeypatch)
    lib = paths["negative"].resolve()
    assert plan["library_mode"] == "msp_only"
    assert plan["lbm"] == {"path": None, "sha256": None}
    assert plan["msp"]["path"] == str(lib)
    assert plan["msp"]["file"] == "lab_neg.msp"
    assert plan["msp"]["source"] == "config_file"
    assert len(plan["msp"]["sha256"]) == 64
    assert plan["removed_declarations"] == [method_file_mod.LBM_KEY]
    assert plan["method"]["overrides"][method_file_mod.LBM_KEY] == ""
    assert plan["method"]["overrides"][method_file_mod.MSP_KEY] == str(lib)

    snapshot = stage_inputs(plan, tmp_path / "pipeline_run")
    effective = (tmp_path / "pipeline_run" / "inputs" / "effective-method.txt").read_text(encoding="ascii")
    assert "Lbm file path: \n" in effective
    assert f"Msp file path: {lib}\n" in effective
    verify_inputs(snapshot)


def test_msp_only_ignores_an_unresolvable_lbm_declaration(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch)
    src = make_source(tmp_path / "raw")
    src["method"].write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\nLbm file path: missing.lbm2\n",
        encoding="ascii", newline="\n")
    request = resolve_request(src["root"], {"library_mode": "msp_only"})
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    assert plan["lbm"]["path"] is None


def test_verify_inputs_detects_a_changed_msp(tmp_path, monkeypatch):
    src, paths, plan = _msp_only_plan(tmp_path, monkeypatch)
    snapshot = stage_inputs(plan, tmp_path / "pipeline_run")
    paths["negative"].write_text("NAME: changed\n", encoding="ascii")
    with pytest.raises(DomainError, match="INPUT_CHANGED") as info:
        verify_inputs(snapshot)
    assert info.value.details["which"] == "msp"


def test_msp_only_without_a_configured_library_stops(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"], {"library_mode": "msp_only"})
    with pytest.raises(DomainError, match="MSP_NOT_CONFIGURED"):
        inspect_inputs(src["root"], request, exe_path=src["exe"])


def test_msp_only_unreadable_library_stops_with_an_envelope(tmp_path, monkeypatch):
    """Review Focus 1: NAS が一時的に読めないときは封筒で止まる（例外で落ちない）。"""
    _allow_fake_exe(monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    src = make_source(tmp_path / "raw")
    request = resolve_request(src["root"], {"library_mode": "msp_only"})
    from metabolomix.pipeline import inputs as inputs_mod
    real = inputs_mod._sha256_file

    def flaky(path):
        if Path(path).suffix == ".msp":
            raise OSError("network name is no longer available")
        return real(path)

    monkeypatch.setattr(inputs_mod, "_sha256_file", flaky)
    with pytest.raises(DomainError, match="MSP_NOT_FOUND") as info:
        inspect_inputs(src["root"], request, exe_path=src["exe"])
    assert "secret_share" not in str(info.value)


def test_msp_only_accepts_a_capitalized_polarity(tmp_path, monkeypatch):
    """Review Focus 3: 極性の値が大文字始まりでも設定キーを引ける。"""
    _allow_fake_exe(monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch)
    src = make_source(tmp_path / "raw")
    src["method"].write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\nLbm file path: fake.lbm2\n",
        encoding="ascii", newline="\n")
    request = resolve_request(src["root"], {"library_mode": "msp_only"})
    plan = inspect_inputs(src["root"], request, exe_path=src["exe"])
    assert plan["msp"]["file"] == "lab_neg.msp"


def test_auto_plan_reports_msp_not_used(tmp_path, monkeypatch):
    _allow_fake_exe(monkeypatch)
    src = make_source(tmp_path / "raw")
    plan = inspect_inputs(src["root"], resolve_request(src["root"]), exe_path=src["exe"])
    assert plan["library_mode"] == "auto"
    assert plan["msp"] == {"path": None, "file": None, "sha256": None, "source": "not_used"}
    assert plan["removed_declarations"] == []
    assert plan["lbm"]["path"] == str(src["lbm"].resolve())
```

- [ ] **Step 2: 失敗するテストを書く（指紋・receipt・レポート）**

`tests/test_pipeline_service.py` の `test_plan_receipt_shows_the_resolved_settings_the_instructions_promise` の
`assert set(resolved) == {"method", "lbm", "polarity"}` を次に置き換える（意図した変更: R4 の MSP を載せる）:

```python
    assert set(resolved) == {"method", "lbm", "polarity", "library_mode", "msp"}
    assert resolved["library_mode"] == "auto"
    assert resolved["msp"] == {"file": None, "sha256": None, "source": "not_used"}
```

同じファイルの末尾に追加:

```python
def test_plan_receipt_shows_the_lab_msp_without_its_location(tmp_path, monkeypatch):
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_only"})
    resolved = receipt["resolved"]
    assert resolved["library_mode"] == "msp_only"
    assert resolved["lbm"] == {"path": None, "sha256": None}
    assert resolved["msp"]["file"] == "lab_neg.msp"
    assert len(resolved["msp"]["sha256"]) == 64
    assert resolved["msp"]["source"] == "config_file"
    assert "secret_share" not in json.dumps(receipt, ensure_ascii=False)


def test_plan_fingerprint_is_unchanged_without_an_msp():
    from metabolomix.core.atomic_io import canonical_hash
    plan = {
        "raw_stat": [{"relative_path": "a.wiff", "size": 1, "mtime_ns": 2}],
        "selected_format": "wiff",
        "method": {"sha256": "m" * 64}, "lbm": {"sha256": "l" * 64},
        "exe": {"sha256": "e" * 64}, "polarity": {"value": "negative"},
    }
    legacy = canonical_hash({
        "raw_stat": plan["raw_stat"], "selected_format": "wiff",
        "method_sha256": "m" * 64, "lbm_sha256": "l" * 64, "exe_sha256": "e" * 64,
        "polarity": "negative"})
    assert service._plan_fingerprint(plan) == legacy
    plan["msp"] = {"path": None, "file": None, "sha256": None, "source": "not_used"}
    assert service._plan_fingerprint(plan) == legacy
    plan["msp"] = {"path": "x", "file": "x.msp", "sha256": "s" * 64, "source": "config_file"}
    assert service._plan_fingerprint(plan) != legacy
```

`tests/test_pipeline_report.py` の末尾に追加:

```python
def test_method_section_shows_the_msp_file_but_not_its_location():
    from metabolomix.pipeline.report import _section_method
    record = {"inputs": {
        "library_mode": "msp_only",
        "method": {"source_path": "m.txt", "sha256": "a" * 64},
        "lbm": {"path": None, "sha256": None},
        "msp": {"path": "//nas/secret_share/lab_neg.msp", "file": "lab_neg.msp",
                "sha256": "b" * 64, "source": "config_file"},
        "exe": {"version": None},
    }}
    text = _section_method(record)
    assert "- library_mode: msp_only" in text
    assert "- msp_file: lab_neg.msp" in text
    assert f"- msp_sha256: {'b' * 64}" in text
    assert "- msp_source: config_file" in text
    assert "secret_share" not in text
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pipeline_inputs.py tests/test_pipeline_service.py tests/test_pipeline_report.py -q -k "msp or resolved_settings or fingerprint or library"`
Expected: FAIL（`KeyError: 'library_mode'`、`resolved` のキー集合の不一致など）

- [ ] **Step 4: 実装する（inputs）**

`metabolomix/pipeline/inputs.py` — `_resolve_lbm_pinned` の直後に追加:

```python
_MSP_NOT_USED = {"path": None, "file": None, "sha256": None, "source": "not_used"}


def _resolve_msp_pinned(request: dict, polarity: str, source_root: Path) -> dict:
    """library_mode="msp_only" の研究室 MSP を解決し、sha256 で固定する。

    エラー文と details に置き場所（ディレクトリ）を載せない（ファイル名だけ。
    core/path_resolvers.LibraryPathError と同じ規約）。絶対パスは入力計画（データ側）にだけ置く。
    """
    override = request.get("msp_file")
    if override and not Path(override).is_absolute():
        override = str(source_root / override)
    polarity = polarity.lower()
    key = user_config.MSP_SETTING_KEYS[polarity]
    try:
        msp_setting = user_config.get_setting(key)
    except user_config.ConfigInvalidError as exc:
        raise DomainError(exc.code, exc.message, exc.details()) from exc
    msp = method_file_mod.resolve_msp(override, polarity, msp_setting)
    if msp.error_code:
        details = ({"msp_file": Path(override).name} if msp.source == "argument"
                   else user_config.describe_missing(key, msp_setting))
        raise DomainError(msp.error_code, msp.message or "", details)
    resolved = Path(msp.path).resolve()
    if not str(resolved).isascii():
        raise DomainError(
            "METHOD_ENCODING_UNSUPPORTED",
            f"研究室 MSP（{resolved.name}）のパスに ASCII 以外の文字があります。MS-DIAL Console は"
            "メソッドファイルを ASCII で読むため、このライブラリは黙って無視されます。",
            {"keys": [method_file_mod.MSP_KEY]})
    try:
        digest = _sha256_file(resolved)
    except OSError as exc:
        raise DomainError(
            "MSP_NOT_FOUND",
            f"研究室 MSP（{resolved.name}）を読めません（ネットワークの一時的な切断の可能性）: "
            f"{type(exc).__name__}",
            {"msp_file": resolved.name}) from exc
    return {"path": str(resolved), "file": resolved.name, "sha256": digest, "source": msp.source}
```

`inspect_inputs` を変える:
- `method_file_mod.resolve_method_references(method_keys, method_path)` の行を次に置き換える:

```python
    library_mode = request.get("library_mode", "auto")
    # 最終選択後にだけ厳格な参照解決を行う（select_methodのグルーピングは弱い版）。
    # msp_only では LBM を使わないので、その宣言の実在で止めない（spec 2026-10-07 §3.2）。
    if library_mode != "msp_only":
        method_file_mod.resolve_method_references(method_keys, method_path)
```

（元のコメント行「# 最終選択後にだけ…」は上に移したので重複させない。）
- `lbm_info = _resolve_lbm_pinned(...)` の行を次に置き換える:

```python
    if library_mode == "msp_only":
        lbm_info = {"path": None, "sha256": None, "source": "disabled"}
        msp_info = _resolve_msp_pinned(request, polarity["value"], source_root)
    else:
        lbm_info = _resolve_lbm_pinned(method_keys, method_path, exe_info["path"], request, source_root)
        msp_info = dict(_MSP_NOT_USED)
```

- `overrides = method_file_mod.relative_path_overrides(method_keys, method_path)` の後の `if lbm_info["path"] and (...)` のブロックを次で囲む:

```python
    removed_declarations: list[str] = []
    if library_mode == "msp_only":
        msp_overrides, removed_declarations = method_file_mod.msp_only_overrides(
            method_keys, msp_info["path"])
        overrides.update(msp_overrides)
    elif lbm_info["path"] and (lbm_info["source"] != "method_file"
                               or method_file_mod.LBM_KEY in overrides):
        # （既存のコメント 3 行をそのまま残す）
        overrides[method_file_mod.LBM_KEY] = lbm_info["path"]
```

- 戻り値の dict の `"lbm": {...},` の直後に追加:

```python
        "library_mode": library_mode,
        "msp": msp_info,
        "removed_declarations": removed_declarations,
```

`verify_inputs` の `for label, info, path_key in (...)` の組に `("msp", snapshot.get("msp"), "path"),` を `("lbm", ...)` の直後に足す。

- [ ] **Step 5: 実装する（service・report・文書）**

`metabolomix/pipeline/service.py` `_plan_fingerprint` の `return canonical_hash({...})` を次に置き換える:

```python
    content = {
        "raw_stat": plan["raw_stat"],
        "selected_format": plan["selected_format"],
        "method_sha256": plan["method"]["sha256"],
        "lbm_sha256": plan["lbm"]["sha256"],
        "exe_sha256": plan["exe"]["sha256"],
        "polarity": plan["polarity"]["value"],
    }
    # 研究室 MSP は msp_only のときだけ入れる。auto の計画の指紋は変えない
    # （受付冪等性と再開の同一性が既存 run で崩れないように。spec 2026-10-07 §3.3）。
    msp_sha256 = (plan.get("msp") or {}).get("sha256")
    if msp_sha256:
        content["msp_sha256"] = msp_sha256
    return canonical_hash(content)
```

`_resolved_settings` の戻り値を次に置き換え、docstring の「載せるのは3項目だけ」を「載せるのは method / lbm / polarity / library_mode / msp だけ（msp は置き場所を載せずファイル名・sha256・出どころ）」に直す:

```python
    msp = inputs.get("msp") or {}
    return {
        "method": {"source_path": method.get("source_path"), "sha256": method.get("sha256")},
        "lbm": {"path": lbm.get("path"), "sha256": lbm.get("sha256")},
        "polarity": {"value": polarity.get("value"), "source": polarity.get("source")},
        "library_mode": inputs.get("library_mode", "auto"),
        "msp": {"file": msp.get("file"), "sha256": msp.get("sha256"),
                "source": msp.get("source", "not_used")},
    }
```

`metabolomix/pipeline/report.py` `_section_method` の `lines = [...]` の `f"- lbm_sha256: ..."` の直後に追加（`msp = inputs.get("msp") or {}` を `exe = ...` の後に足す）:

```python
        f"- library_mode: {inputs.get('library_mode') or 'auto'}",
        f"- msp_file: {msp.get('file') or '(not used)'}",
        f"- msp_sha256: {msp.get('sha256') or '(not used)'}",
        f"- msp_source: {msp.get('source') or 'not_used'}",
```

`metabolomix/tools/pipeline_tools.py` `pipeline_plan` の docstring:
- `request:` の節の項目列挙に `library_mode` / `msp_file` を足し、次の 2 文を足す:
  「`library_mode="msp_only"` は LBM を使わず研究室の参照ライブラリ（.msp）だけで同定する（`msp_file` か極性の設定 `[library] msp_positive` / `msp_negative`。メソッドの Text DB・注釈器設定表の宣言は空にする）。`lbm_file` との同時指定、`auto` での `msp_file` は `PIPELINE_REQUEST_INVALID`。」
- receipt の `resolved` の説明を「`method`・`lbm`・`polarity`・`library_mode`・`msp`（ファイル名・sha256・出どころ。置き場所は載せない）」に直す。

`USAGE.md:199`（`pipeline_plan` の行）の「receiptの `resolved` が返す method/lbm/polarity」を「receiptの `resolved` が返す method/lbm/polarity/library_mode/msp」に直し、行末に「要求の `library_mode="msp_only"` で LBM を使わず研究室 MSP だけの解析になる(`console_plan` と同じ規則。MSP の sha256 を計画時と実行時に確かめ、変わっていれば `INPUT_CHANGED`)。」を足す。

`docs/workflow/pipeline.md`: 入力計画（method / LBM / 実行体を決める段）を説明している節の末尾に、次の段落を足す（関数名を書くなら実在するものだけ。`tests/test_workflow_docs.py` が検証する）:

```markdown
`library_mode="msp_only"` の要求では LBM を解決せず（`Lbm file path` は空）、研究室の参照ライブラリ（.msp）を
`msp_file` か極性の設定（`[library] msp_positive` / `msp_negative`）から決めて実効メソッドの `Msp file path` に書く。
メソッドが宣言する Text DB と注釈器設定表は空にする。MSP の sha256 は入力計画に固定し、実行時に `verify_inputs` が
再検査する（変わっていれば `INPUT_CHANGED`、`which: "msp"`）。receipt にはファイル名・sha256・出どころだけを載せ、
置き場所は載せない。
```

- [ ] **Step 6: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_pipeline_inputs.py tests/test_pipeline_service.py tests/test_pipeline_report.py tests/test_pipeline_recovery.py tests/test_workflow_docs.py tests/test_readme_links.py -q`
Expected: PASS（全件）

- [ ] **Step 7: コミットする（バックグラウンド）**

```bash
cd /c/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-library-mode-msp-only
git add metabolomix/pipeline/inputs.py metabolomix/pipeline/service.py metabolomix/pipeline/report.py metabolomix/tools/pipeline_tools.py USAGE.md docs/workflow/pipeline.md tests/test_pipeline_inputs.py tests/test_pipeline_service.py tests/test_pipeline_report.py
git commit -m "feat(pipeline): pipeline v1 に library_mode=msp_only を足す（MSP を sha256 で固定し実行時に再検査）

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$TMP/commit_task4.log" 2>&1; echo "exit=$?"; grep -E "passed|failed|FAILED" "$TMP/commit_task4.log" | tail -5
```

Expected: `exit=0`、全テスト通過。

---

## 実装後（全 Task の完了後、最終レビューの後）

1. main ツリーの `docs/HISTRY.md` に日付見出し `## 2026-10-07 library_mode="msp_only"（reanalysis-study Task 9）` で、何を足したか・確かめた事実・意図して変えた既存テスト（`test_defaults_are_fully_populated` に 2 行、`resolved` のキー集合）を追記する。`docs/task.md` から「LBM 無効化・研究室 MSP 指定（Task 9）」と `pipeline/Console の Msp file path を MSDIAL_MSP_* から自動で埋める件` の項目を消す（並列のエージェントがいないことを確かめてから）。
2. main ツリーで `git merge --no-ff feat/library-mode-msp-only -m "Merge feat/library-mode-msp-only: Console と pipeline v1 に library_mode=msp_only（LBM 不使用・研究室 MSP のみ）を足す"`（pre-commit 相当の全テストはマージ後に 1 回走らせて確かめる）。push はユーザーに確認する。
3. reanalysis-study の手順書 工程 5 を確定する（reanalysis-study plan Task 9 Step 4〜5。`library_mode="msp_only"` を `pipeline_plan` / `pipeline_run` の要求に書き、版を v3 に上げる）。
