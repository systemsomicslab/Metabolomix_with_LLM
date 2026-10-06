# 外部資産のパスを設定ファイルから読む Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** MS-DIAL Console の実行体・脂質ライブラリ・研究室の参照ライブラリの場所を、環境変数に加えてリポジトリ直下の `lipidmix.local.toml` から読めるようにする。

**Architecture:** stdlib だけの leaf モジュール `lipidmix/core/user_config.py` が「環境変数 → 設定ファイル」の優先順と TOML の解析・検証・エラー文の組み立てを一手に持つ。既存のリゾルバ（`console/runner.get_exe_path`、`console/method_file.resolve_lbm`、`core/path_resolvers` の研究室ライブラリ解決）は自分で `os.environ` を読むのをやめ、このモジュールから値を引く。エラーコードは据え置き、`details` を足す。

**Tech Stack:** Python 3.13+（開発は 3.14）、stdlib `tomllib`、pytest。

**Spec:** `docs/superpowers/specs/2026-10-06-user-config-file-design.md`（実装者は計画と併せて読む）

## Global Constraints

- 設定ファイル: `LIPIDMIX_CONFIG` が指すファイル → `<repo>/lipidmix.local.toml` → 無ければ設定なし。1 ファイルだけ読み、マージしない。
- 雛形は `lipidmix.example.toml`（追跡対象）。`lipidmix.local.toml` は `.gitignore` に入れる。
- キーは 4 つだけ: `[msdial] exe`（`MSDIAL_EXE`）、`[msdial] lbm`（`MSDIAL_LBM`）、`[library] msp_positive`（`MSDIAL_MSP_POS`）、`[library] msp_negative`（`MSDIAL_MSP_NEG`）。
- 値を引く順: 環境変数 → 設定ファイル → 各リゾルバの既存の推定。空文字・空白だけの環境変数は未設定。環境変数があるキーでは設定ファイルを読まない。
- 環境変数の値は前後の空白を除いてそのまま返す（絶対パス化しない。既存テストが `"fake.exe"` や `"C:/MsDial/..."` を置いている）。設定ファイルの値は `~` を展開し、相対ならそのファイルのフォルダ基準で絶対パスにする。
- 設定ファイルは `utf-8-sig` で読む（BOM 可）。呼ばれるたびに読み、`(パス, st_mtime_ns, st_size)` で解析結果をキャッシュする（再起動不要）。
- `lipidmix/core/user_config.py` は stdlib だけを import する（`mcp_core` / `session_state` / `lipidmix.tools.*` 不可）。
- エラーコードは変えない: `MSDIAL_EXE_NOT_FOUND` / `LBM_NOT_FOUND` / `MSP_ENV_NOT_FOUND`。新設は `CONFIG_INVALID` だけ。
- 研究室ライブラリ（`msp_*`）の値は、エラー文にも `details` にもディレクトリを出さない（ファイル名とキー名だけ）。設定ファイル自身のパスは出してよい。
- 設定の不備でサーバの起動を止めない。止まるのはその設定を要するツールだけ。
- `required_tools` は付けない。
- README.md / CLAUDE.md に数量表現を書かない（`tests/test_readme_links.py`）。
- コミットは pre-commit で全テストが走る（数分）。`git commit` はバックグラウンドで実行し、出力をファイルに落として `FAILED` を確認する。実行中に作業ツリーを編集しない。
- 作業場所は worktree `C:\Users\yuu18\Metabolomix_with_LLM\.worktrees\feat-user-config-file`（branch `feat/user-config-file`）。Python は `C:/Python314/python.exe`。

## Review Focus

1. **メモ帳で保存した設定ファイル（UTF-8 BOM 付き）**: 普通に読めること。`tomllib` は BOM を `Invalid statement` にするので、`utf-8-sig` で読まないと外部の利用者が最初の一歩で詰まる。→ Task 1 `test_utf8_bom_is_accepted`。
2. **二重引用符で書いた Windows パスのうち、構文として通ってしまうもの**（`"C:\new\tool.exe"`）: `\n` `\t` が改行・タブになって黙って壊れたパスになる。`CONFIG_INVALID` で単一引用符を案内すること。→ Task 1 `test_escape_that_parses_into_control_chars_is_config_invalid`。
3. **設定ファイルが壊れているが、環境変数で値を渡している既存環境**: 今まで通り動くこと（設定ファイルを読みに行かない）。→ Task 1 `test_env_var_bypasses_a_broken_config_file`。
4. **研究室ライブラリの場所の漏洩**: 設定ファイル由来の `msp_*` が指す先が無いとき、`library_load` の戻り値にディレクトリが出ないこと。→ Task 4 `test_library_load_missing_config_msp_does_not_leak_the_directory`。
5. **exe は環境変数・LBM は設定ファイル、のように出どころが混ざる構成**: 各値がそれぞれの出どころで解決され、LBM の `source` が `config_file` になること。→ Task 3 `test_template_uses_lbm_from_the_config_file_while_exe_comes_from_env`。

---

### Task 1: 設定モジュール `user_config` と、テストの隔離

**Files:**
- Create: `lipidmix/core/user_config.py`
- Create: `tests/test_user_config.py`
- Modify: `tests/conftest.py`（全体を置き換え）
- Modify: `.gitignore`（3 行目 `.env` の直後に 1 行）
- Modify: `tests/test_gitignore_library.py`（`_MUST_BE_IGNORED`）

**Interfaces:**
- Consumes: なし
- Produces（後続タスクが使う。名前と型を変えない）:
  - `user_config.CONFIG_ENV = "LIPIDMIX_CONFIG"`、`CONFIG_FILENAME = "lipidmix.local.toml"`、`EXAMPLE_FILENAME = "lipidmix.example.toml"`、`REPO_ROOT: Path`
  - `user_config.SETTINGS: dict[str, str]`（`"msdial.exe" → "MSDIAL_EXE"` ほか 4 件）
  - `@dataclass(frozen=True) class Setting: key: str; value: str; source: str  # "env" | "config_file"; env_var: str; config_file: str | None`
  - `class ConfigInvalidError(Exception)`: 属性 `code = "CONFIG_INVALID"`、`message: str`、`config_file: str`、`line: int | None`、`column: int | None`、メソッド `details() -> dict`
  - `config_file_path() -> Path`
  - `get_setting(key: str) -> Setting | None`（設定ファイルが読めなければ `ConfigInvalidError`）
  - `setting_label(setting: Setting) -> str`
  - `missing_hint(key: str, what: str) -> str`（送出しない）
  - `describe_missing(key: str, setting: Setting | None = None) -> dict`（送出しない）

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_user_config.py`:

```python
"""外部資産の場所の設定（`lipidmix.core.user_config`）。spec 2026-10-06。

環境変数 → 設定ファイル（`LIPIDMIX_CONFIG` → `<repo>/lipidmix.local.toml`）の順に引く。
conftest が環境変数を消し `LIPIDMIX_CONFIG` を存在しないパスに向けているので、
各テストは必要なものだけを置く。
"""
from __future__ import annotations

import ast
import json
import os
import sys
from pathlib import Path

import pytest

from lipidmix.core import user_config
from lipidmix.core.user_config import ConfigInvalidError

BS = "\\"


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    """tmp の設定ファイル（まだ書かない）を LIPIDMIX_CONFIG で指す。"""
    path = tmp_path / "conf" / "lipidmix.local.toml"
    path.parent.mkdir()
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(path))
    return path


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


# ---------- 探す順 ----------

def test_config_file_path_follows_lipidmix_config(cfg):
    assert user_config.config_file_path() == cfg


def test_config_file_path_defaults_to_the_repository_root(monkeypatch):
    from lipidmix.core import mcp_core
    monkeypatch.delenv("LIPIDMIX_CONFIG", raising=False)
    assert user_config.REPO_ROOT == mcp_core.BASE_DIR
    assert user_config.config_file_path() == mcp_core.BASE_DIR / "lipidmix.local.toml"


def test_absent_config_file_means_no_setting(cfg):
    assert user_config.get_setting("msdial.exe") is None


# ---------- 引く順と値の形 ----------

def test_reads_a_value_from_the_config_file(cfg, tmp_path):
    exe = tmp_path / "MSDIALCUI.exe"
    _write(cfg, f"[msdial]\nexe = '{exe}'\n")
    setting = user_config.get_setting("msdial.exe")
    assert setting == user_config.Setting(
        key="msdial.exe", value=str(exe), source="config_file",
        env_var="MSDIAL_EXE", config_file=str(cfg))


def test_env_var_wins_over_the_config_file(cfg, tmp_path, monkeypatch):
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'from_file.exe'}'\n")
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    setting = user_config.get_setting("msdial.exe")
    assert (setting.value, setting.source, setting.config_file) == ("fake.exe", "env", None)


def test_blank_env_var_counts_as_unset(cfg, tmp_path, monkeypatch):
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'from_file.exe'}'\n")
    monkeypatch.setenv("MSDIAL_EXE", "   ")
    assert user_config.get_setting("msdial.exe").source == "config_file"


def test_env_value_is_returned_verbatim(monkeypatch):
    monkeypatch.setenv("MSDIAL_EXE", " C:/MsDial/MsdialConsoleApp.exe ")
    assert user_config.get_setting("msdial.exe").value == "C:/MsDial/MsdialConsoleApp.exe"


def test_relative_path_resolves_against_the_config_folder(cfg):
    _write(cfg, "[library]\nmsp_positive = 'libs/pos.msp'\n")
    setting = user_config.get_setting("library.msp_positive")
    assert setting.value == str(cfg.parent / "libs" / "pos.msp")


def test_tilde_is_expanded(cfg, tmp_path, monkeypatch):
    home = tmp_path / "home"
    monkeypatch.setenv("USERPROFILE", str(home))
    monkeypatch.setenv("HOME", str(home))
    _write(cfg, "[msdial]\nexe = '~/MSDIALCUI.exe'\n")
    assert user_config.get_setting("msdial.exe").value == str(home / "MSDIALCUI.exe")


def test_empty_value_counts_as_unset(cfg):
    _write(cfg, "[msdial]\nexe = '  '\n")
    assert user_config.get_setting("msdial.exe") is None


def test_rewritten_file_takes_effect_without_restart(cfg, tmp_path):
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'a.exe'}'\n")
    assert user_config.get_setting("msdial.exe").value == str(tmp_path / "a.exe")
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'bb.exe'}'\n")
    stat = cfg.stat()
    os.utime(cfg, ns=(stat.st_atime_ns, stat.st_mtime_ns + 2_000_000_000))
    assert user_config.get_setting("msdial.exe").value == str(tmp_path / "bb.exe")


def test_utf8_bom_is_accepted(cfg, tmp_path):
    """メモ帳は BOM を付けることがある。tomllib は BOM を構文エラーにする。"""
    cfg.write_text(f"\ufeff[msdial]\nexe = '{tmp_path / 'a.exe'}'\n", encoding="utf-8")
    assert user_config.get_setting("msdial.exe").value == str(tmp_path / "a.exe")


# ---------- 読めない設定ファイル ----------

def test_double_quoted_windows_path_is_config_invalid_with_a_hint(cfg):
    _write(cfg, '[msdial]\nexe = "C:' + BS + 'MS-DIAL' + BS + 'MSDIALCUI.exe"\n')
    with pytest.raises(ConfigInvalidError) as exc:
        user_config.get_setting("msdial.exe")
    assert exc.value.code == "CONFIG_INVALID"
    assert "単一引用符" in exc.value.message
    assert exc.value.details()["config_file"] == str(cfg)
    if sys.version_info >= (3, 14):  # TOMLDecodeError.lineno は 3.14 から
        assert exc.value.line == 2


def test_escape_that_parses_into_control_chars_is_config_invalid(cfg):
    """`\\n` `\\t` は TOML の正しいエスケープなので構文エラーにならず、黙って壊れる。"""
    _write(cfg, '[msdial]\nexe = "C:' + BS + 'new' + BS + 'tool.exe"\n')
    with pytest.raises(ConfigInvalidError) as exc:
        user_config.get_setting("msdial.exe")
    assert "単一引用符" in exc.value.message


def test_non_string_value_is_config_invalid(cfg):
    _write(cfg, "[msdial]\nexe = 3\n")
    with pytest.raises(ConfigInvalidError) as exc:
        user_config.get_setting("msdial.exe")
    assert "文字列" in exc.value.message


def test_env_var_bypasses_a_broken_config_file(cfg, monkeypatch):
    _write(cfg, "[msdial\n")
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    assert user_config.get_setting("msdial.exe").value == "fake.exe"


# ---------- エラー文の材料 ----------

def test_describe_missing_without_a_config_file(cfg):
    assert user_config.describe_missing("msdial.exe") == {
        "setting": "msdial.exe", "env_var": "MSDIAL_EXE", "source": None,
        "config_file": str(cfg), "config_file_exists": False,
        "example": "lipidmix.example.toml"}


def test_unknown_keys_are_reported(cfg):
    _write(cfg, "[msdail]\nexe = 'a.exe'\n[library]\nmsp_postive = 'b.msp'\n")
    assert user_config.get_setting("msdial.exe") is None
    details = user_config.describe_missing("msdial.exe")
    assert details["unknown_keys"] == ["library.msp_postive", "msdail.exe"]
    assert "msdail.exe" in user_config.missing_hint("msdial.exe", "MS-DIAL Console の実行体")


def test_describe_missing_never_contains_the_value(cfg, tmp_path):
    secret = tmp_path / "secret-lab-share" / "pos.msp"
    _write(cfg, f"[library]\nmsp_positive = '{secret}'\n")
    setting = user_config.get_setting("library.msp_positive")
    details = user_config.describe_missing("library.msp_positive", setting)
    assert details["source"] == "config_file"
    assert "secret-lab-share" not in json.dumps(details, ensure_ascii=False)


def test_describe_missing_reports_a_broken_file_without_raising(cfg):
    _write(cfg, "[msdial\n")
    assert "TOML" in user_config.describe_missing("msdial.exe")["config_invalid"]


def test_missing_hint_when_the_config_file_is_absent(cfg):
    hint = user_config.missing_hint("msdial.exe", "MS-DIAL Console の実行体")
    assert "lipidmix.example.toml" in hint
    assert "[msdial] exe" in hint
    assert "MSDIAL_EXE" in hint


def test_missing_hint_when_the_config_file_exists(cfg):
    _write(cfg, "[library]\n")
    hint = user_config.missing_hint("msdial.exe", "MS-DIAL Console の実行体")
    assert hint.startswith("lipidmix.local.toml の [msdial] exe に")


def test_setting_label_names_the_source(cfg, tmp_path, monkeypatch):
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'a.exe'}'\n")
    assert user_config.setting_label(user_config.get_setting("msdial.exe")) == \
        "lipidmix.local.toml の [msdial] exe"
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    assert user_config.setting_label(user_config.get_setting("msdial.exe")) == \
        "環境変数 MSDIAL_EXE"


# ---------- leaf ----------

def test_module_imports_only_the_standard_library():
    """pipeline の独立 worker も import する。重い依存・グローバル session を持ち込まない。"""
    source = Path(user_config.__file__).read_text(encoding="utf-8")
    roots = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0:
            roots.add(node.module.split(".")[0])
    assert roots - {"__future__"} <= set(sys.stdlib_module_names)
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_user_config.py -q`
Expected: 収集時に `ModuleNotFoundError: No module named 'lipidmix.core.user_config'`

- [ ] **Step 3: `lipidmix/core/user_config.py` を書く**

```python
"""外部資産（MS-DIAL Console・脂質ライブラリ・研究室の参照ライブラリ）の場所の設定。

spec: docs/superpowers/specs/2026-10-06-user-config-file-design.md

値を引く順は 環境変数 → 設定ファイル。設定ファイルは `LIPIDMIX_CONFIG` が指す
ファイル、無ければ `<repo>/lipidmix.local.toml`（雛形 `lipidmix.example.toml`）。
環境変数が設定されているキーでは設定ファイルを読まない——設定ファイルが壊れて
いても、環境変数で渡している既存環境は今まで通り動く。

呼ばれるたびに読む（パス・更新時刻・大きさが同じなら前回の解析結果を使う）ので、
設定ファイルの書き換えはサーバを再起動しなくても次の呼び出しから効く。

**stdlib だけの leaf**。pipeline の独立 worker も import するので、`mcp_core` /
`session_state` / `lipidmix.tools.*` を import しない（tests/test_user_config.py が縛る）。
"""
from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path

CONFIG_ENV = "LIPIDMIX_CONFIG"
CONFIG_FILENAME = "lipidmix.local.toml"
EXAMPLE_FILENAME = "lipidmix.example.toml"
#: core/data_config.py の DEFAULT_DATA_DIR と同じ起点（mcp_core.BASE_DIR との一致は
#: tests/test_user_config.py が縛る。mcp_core は leaf ではないのでここから参照しない）。
REPO_ROOT = Path(__file__).resolve().parents[2]

#: 設定キー（`<節>.<名前>`）→ 同じ値を指す環境変数。環境変数名の正準はここ。
SETTINGS: dict[str, str] = {
    "msdial.exe": "MSDIAL_EXE",
    "msdial.lbm": "MSDIAL_LBM",
    "library.msp_positive": "MSDIAL_MSP_POS",
    "library.msp_negative": "MSDIAL_MSP_NEG",
}

# 二重引用符の中の `\n` `\t` `\b` `\f` `\r` は TOML の正しいエスケープなので構文エラーに
# ならない。`"C:\new\tool.exe"` は黙って改行とタブ入りの値になる。値の側で拾う。
_CONTROL_CHARS = frozenset("\n\t\b\f\r")
_QUOTE_HINT = ("Windows のパスは単一引用符で囲んでください"
               "（例: exe = 'C:\\MS-DIAL\\MSDIALCUI.exe'）。")


@dataclass(frozen=True)
class Setting:
    """1 つの設定値と、その出どころ。"""

    key: str
    #: 環境変数の値は前後の空白を除いてそのまま（既存の利用者とテストが相対名や
    #: `/` 区切りを置いているので形を変えない）。設定ファイルの値は絶対パス。
    value: str
    source: str  # "env" | "config_file"
    env_var: str
    config_file: str | None  # source == "config_file" のときの設定ファイル


class ConfigInvalidError(Exception):
    """設定ファイルが読めない（構文・型・壊れたパス）。メッセージは利用者向け。"""

    code = "CONFIG_INVALID"

    def __init__(self, message: str, *, config_file: str,
                 line: int | None = None, column: int | None = None):
        super().__init__(message)
        self.message = message
        self.config_file = config_file
        self.line = line
        self.column = column

    def details(self) -> dict:
        return {"config_file": self.config_file, "line": self.line, "column": self.column}


@dataclass(frozen=True)
class _Parsed:
    values: dict[str, str]
    unknown_keys: tuple[str, ...]


_cache: dict[tuple[str, int, int], _Parsed] = {}


def config_file_path() -> Path:
    """読むべき設定ファイルのパス（存在しなくても返す）。"""
    override = (os.environ.get(CONFIG_ENV) or "").strip()
    if override:
        return Path(override).expanduser()
    return REPO_ROOT / CONFIG_FILENAME


def _resolve(value: str, base: Path) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return str(path)


def _parse(path: Path) -> _Parsed | None:
    """設定ファイルを解析する。無ければ None、読めなければ ConfigInvalidError。"""
    if not path.is_file():
        return None
    stat = path.stat()
    cache_key = (str(path), stat.st_mtime_ns, stat.st_size)
    cached = _cache.get(cache_key)
    if cached is not None:
        return cached

    try:
        text = path.read_text(encoding="utf-8-sig")  # メモ帳の BOM を許す
    except (OSError, UnicodeDecodeError) as exc:
        raise ConfigInvalidError(f"{path.name} を UTF-8 として読めません: {exc}",
                                 config_file=str(path)) from exc
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        message = f"{path.name} の TOML 構文エラー: {exc}"
        if "\\" in text:
            message += f" {_QUOTE_HINT}"
        raise ConfigInvalidError(message, config_file=str(path),
                                 line=getattr(exc, "lineno", None),
                                 column=getattr(exc, "colno", None)) from exc

    values: dict[str, str] = {}
    unknown: list[str] = []
    for section, body in data.items():
        if not isinstance(body, dict):
            unknown.append(section)
            continue
        for name, value in body.items():
            key = f"{section}.{name}"
            if key not in SETTINGS:
                unknown.append(key)
                continue
            if not isinstance(value, str):
                raise ConfigInvalidError(
                    f"{path.name} の [{section}] {name} は文字列（パス）で書いてください"
                    f"（受け取った型: {type(value).__name__}）。{_QUOTE_HINT}",
                    config_file=str(path))
            if _CONTROL_CHARS & set(value):
                raise ConfigInvalidError(
                    f"{path.name} の [{section}] {name} に改行やタブが入っています。"
                    f"二重引用符の中の \\n や \\t がエスケープとして読まれた可能性があります。"
                    f"{_QUOTE_HINT}",
                    config_file=str(path))
            value = value.strip()
            if value:
                values[key] = _resolve(value, path.parent)

    parsed = _Parsed(values=values, unknown_keys=tuple(sorted(unknown)))
    _cache.clear()
    _cache[cache_key] = parsed
    return parsed


def get_setting(key: str) -> Setting | None:
    """`key` の値を 環境変数 → 設定ファイル の順に引く。どちらにも無ければ None。

    設定ファイルが読めなければ ConfigInvalidError（環境変数があれば読まないので送出しない）。
    """
    env_var = SETTINGS[key]
    from_env = (os.environ.get(env_var) or "").strip()
    if from_env:
        return Setting(key=key, value=from_env, source="env", env_var=env_var, config_file=None)
    path = config_file_path()
    parsed = _parse(path)
    if parsed is None or key not in parsed.values:
        return None
    return Setting(key=key, value=parsed.values[key], source="config_file",
                   env_var=env_var, config_file=str(path))


def _section_and_name(key: str) -> tuple[str, str]:
    section, name = key.split(".", 1)
    return section, name


def setting_label(setting: Setting) -> str:
    """エラー文で値の出どころを言う句。"""
    if setting.source == "env":
        return f"環境変数 {setting.env_var}"
    section, name = _section_and_name(setting.key)
    return f"{Path(setting.config_file).name} の [{section}] {name}"


def _unknown_keys_quietly(path: Path) -> tuple[tuple[str, ...], str | None]:
    """(未知のキー, 読めなかった理由)。エラー文の材料なので送出しない。"""
    try:
        parsed = _parse(path)
    except ConfigInvalidError as exc:
        return (), exc.message
    return (parsed.unknown_keys if parsed else ()), None


def missing_hint(key: str, what: str) -> str:
    """未設定のとき利用者に伝える文。`what` は「MS-DIAL Console の実行体」など。"""
    section, name = _section_and_name(key)
    path = config_file_path()
    if path.is_file():
        head = f"{path.name} の [{section}] {name} に{what}のパスを書いてください"
    else:
        head = (f"{EXAMPLE_FILENAME} を {path} として複製し、"
                f"[{section}] {name} に{what}のパスを書いてください")
    message = f"{head}（環境変数 {SETTINGS[key]} でも指定できます）。"
    unknown, _ = _unknown_keys_quietly(path)
    if unknown:
        message += f" 認識できないキーがあります: {', '.join(unknown)}。"
    return message


def describe_missing(key: str, setting: Setting | None = None) -> dict:
    """未設定・指す先が無いエラーの `details` に足す内容。値そのものは載せない。"""
    path = config_file_path()
    details: dict = {
        "setting": key,
        "env_var": SETTINGS[key],
        "source": setting.source if setting else None,
        "config_file": str(path),
        "config_file_exists": path.is_file(),
        "example": EXAMPLE_FILENAME,
    }
    unknown, invalid = _unknown_keys_quietly(path)
    if unknown:
        details["unknown_keys"] = list(unknown)
    if invalid:
        details["config_invalid"] = invalid
    return details
```

- [ ] **Step 4: テストの隔離を広げる**

`tests/conftest.py` を次の内容に置き換える:

```python
"""全テスト共通の隔離。

外部資産の場所（`MSDIAL_EXE` / `MSDIAL_LBM` / `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）は、
利用者が OS や `.mcp.json` に常設するか、リポジトリ直下の `lipidmix.local.toml` に書く
前提の設定なので、テストを走らせる機械にも置かれている。残したままだと
`resolve_library_path()` がテストの tmp ではなく実ライブラリを掴み、結果が機械ごとに
変わる（しかも実ライブラリは外部流出禁止の資産）。環境変数は消し、設定ファイルは
存在しないパスへ向ける。使うテストは monkeypatch.setenv で明示的に置く
（設定ファイルは tmp に書いて `LIPIDMIX_CONFIG` で指す）。
"""
from pathlib import Path

import pytest

from lipidmix.core.user_config import CONFIG_ENV, SETTINGS

#: 作られることの無いパス。autouse で tmp_path を要求すると全テストに tmp ディレクトリが
#: できるので、存在しない固定パスで済ませる。
_NO_CONFIG = Path(__file__).resolve().parent / "_no_such_dir" / "lipidmix.local.toml"


@pytest.fixture(autouse=True)
def _isolate_external_asset_settings(monkeypatch):
    for name in SETTINGS.values():
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(CONFIG_ENV, str(_NO_CONFIG))
```

- [ ] **Step 5: `lipidmix.local.toml` を git の網に入れる**

`.gitignore` の 3 行目 `.env` の直後に追加:

```
# 外部資産の場所の設定（研究室の参照ライブラリのパスが入る）。雛形は lipidmix.example.toml
lipidmix.local.toml
```

`tests/test_gitignore_library.py` の `_MUST_BE_IGNORED` の末尾（`f"cache/{store.DIGEST_INDEX_NAME}",` の次）に追加:

```python
    # 設定ファイルには研究室ライブラリのパスが入る（spec 2026-10-06）
    "lipidmix.local.toml",
```

- [ ] **Step 6: 通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_user_config.py tests/test_gitignore_library.py -q`
Expected: 全て PASS

- [ ] **Step 7: 全テストで隔離の副作用を確認する**

Run: `C:/Python314/python.exe -m pytest tests -q -x -p no:cacheprovider 2>&1 | tail -5`
Expected: 全て PASS。もし落ちるテストがあれば、それは手元の `MSDIAL_EXE` / `MSDIAL_LBM` に暗黙に頼っていたテスト。そのテストで `monkeypatch.setenv` を明示して直す（隔離を緩めない）。

- [ ] **Step 8: コミット（バックグラウンド）**

```bash
git add lipidmix/core/user_config.py tests/test_user_config.py tests/conftest.py .gitignore tests/test_gitignore_library.py
git commit -m "feat(config): 外部資産の場所を設定ファイルから引く user_config を足す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$SCRATCH/commit-t1.log" 2>&1; echo exit=$?; grep -E "passed|failed|FAILED" "$SCRATCH/commit-t1.log" | tail -3
```

（`$SCRATCH` はセッションの scratchpad。バックグラウンド実行にし、完了通知を待つ。）

---

### Task 2: Console の実行体（`[msdial] exe`）

**Files:**
- Modify: `lipidmix/console/runner.py`（モジュール docstring 冒頭 3〜4 行目、`MsdialExeNotFoundError`、`get_exe_path`）
- Modify: `lipidmix/tools/console_tools.py`（`console_plan` の exe 取得、`console_run` の exe 取得、`console_method_template` の exe 取得、`_console_run_result` の `launch_failed` 分岐、`_msdial_exe_setup_help`、新規ヘルパ `_exe_error` / `_configured_exe`）
- Modify: `lipidmix/pipeline/service.py`（`_resolve_exe_path`）
- Create: `tests/test_console_config.py`

**Interfaces:**
- Consumes: Task 1 の `user_config.get_setting` / `missing_hint` / `describe_missing` / `config_file_path` / `EXAMPLE_FILENAME` / `ConfigInvalidError`
- Produces:
  - `MsdialExeNotFoundError(message: str, *, code: str = "MSDIAL_EXE_NOT_FOUND", details: dict | None = None)`。`EnvironmentError` のサブクラスのまま。属性 `code: str`、`details: dict`。設定ファイルが読めないときも `code="CONFIG_INVALID"` でこの型を送出する（既存の `except EnvironmentError` をすべて生かすため）。
  - `console_tools._exe_error(exc: EnvironmentError) -> str`（封筒）
  - `console_tools._configured_exe() -> str`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_console_config.py`:

```python
"""Console 経路が設定ファイルから実行体・LBM を引く（spec 2026-10-06 §6・§7）。"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from lipidmix.console.runner import MsdialExeNotFoundError, get_exe_path

BS = "\\"


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    path = tmp_path / "conf" / "lipidmix.local.toml"
    path.parent.mkdir()
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(path))
    return path


def _method(tmp_path: Path) -> Path:
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    return method


def _plan(tmp_path: Path) -> dict:
    from lipidmix.tools.console_tools import console_plan
    return json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(_method(tmp_path)),
                                   polarity="negative", measure="peak_height"))


# ---------- get_exe_path ----------

def test_get_exe_path_reads_the_config_file(cfg, tmp_path):
    exe = tmp_path / "MSDIALCUI.exe"
    cfg.write_text(f"[msdial]\nexe = '{exe}'\n", encoding="utf-8")
    assert get_exe_path() == str(exe)


def test_unset_exe_names_the_config_file_and_the_key(cfg):
    with pytest.raises(MsdialExeNotFoundError) as exc:
        get_exe_path()
    assert exc.value.code == "MSDIAL_EXE_NOT_FOUND"
    assert "lipidmix.example.toml" in str(exc.value)
    assert "[msdial] exe" in str(exc.value)
    assert exc.value.details["config_file"] == str(cfg)
    assert exc.value.details["config_file_exists"] is False


def test_broken_config_file_is_config_invalid(cfg):
    cfg.write_text("[msdial\n", encoding="utf-8")
    with pytest.raises(MsdialExeNotFoundError) as exc:
        get_exe_path()
    assert exc.value.code == "CONFIG_INVALID"
    assert exc.value.details["config_file"] == str(cfg)


# ---------- console_plan の封筒 ----------

def test_console_plan_envelope_carries_the_config_details(tmp_path, cfg):
    err = _plan(tmp_path)["error"]
    assert err["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert err["details"]["human_action_required"] is True
    assert err["details"]["setting"] == "msdial.exe"
    assert err["details"]["restart_required"] is False
    assert any("[msdial] exe" in step for step in err["details"]["how_to_set"])


def test_console_plan_reports_config_invalid_with_the_quote_hint(tmp_path, cfg):
    cfg.write_text('[msdial]\nexe = "C:' + BS + 'MS-DIAL' + BS + 'x.exe"\n', encoding="utf-8")
    err = _plan(tmp_path)["error"]
    assert err["code"] == "CONFIG_INVALID"
    assert "単一引用符" in err["message"]


# ---------- console_method_template ----------

def _neg_param(path: Path) -> Path:
    path.write_text(
        "Ion mode: Negative\nTarget omics: Lipidomics\nLbm file path: \n"
        "Searched adduct ions: [M-H]-\n", encoding="ascii")
    return path


def test_template_uses_the_exe_from_the_config_file_to_find_the_lbm(tmp_path, cfg):
    """exe フォルダの *.lbm2 を 1 件だけ自動採用する GUI 流の推定が、設定ファイルの exe でも効く。"""
    from lipidmix.tools.console_tools import console_method_template
    app = tmp_path / "app"
    app.mkdir()
    (app / "lib.lbm2").touch()
    cfg.write_text(f"[msdial]\nexe = '{app / 'MSDIALCUI.exe'}'\n", encoding="utf-8")
    parsed = json.loads(console_method_template(
        out_path=str(tmp_path / "param_POS.txt"), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["status"] == "written"
    assert parsed["lbm"]["source"] == "exe_dir"


def test_template_reports_config_invalid(tmp_path, cfg):
    from lipidmix.tools.console_tools import console_method_template
    cfg.write_text("[msdial\n", encoding="utf-8")
    parsed = json.loads(console_method_template(
        out_path=str(tmp_path / "param_POS.txt"), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["error"]["code"] == "CONFIG_INVALID"


# ---------- pipeline の受付 ----------

def test_pipeline_intake_resolves_the_exe_from_the_config_file(cfg, tmp_path):
    from lipidmix.pipeline import service
    exe = tmp_path / "MSDIALCUI.exe"
    cfg.write_text(f"[msdial]\nexe = '{exe}'\n", encoding="utf-8")
    assert service._resolve_exe_path() == exe


def test_pipeline_intake_wraps_config_invalid(cfg):
    from lipidmix.core.atomic_io import DomainError
    from lipidmix.pipeline import service
    cfg.write_text("[msdial\n", encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        service._resolve_exe_path()
    assert exc.value.code == "CONFIG_INVALID"


def test_pipeline_intake_keeps_the_not_found_code_with_details(cfg):
    from lipidmix.core.atomic_io import DomainError
    from lipidmix.pipeline import service
    with pytest.raises(DomainError) as exc:
        service._resolve_exe_path()
    assert exc.value.code == "MSDIAL_EXE_NOT_FOUND"
    assert exc.value.details["setting"] == "msdial.exe"
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_console_config.py -q`
Expected: FAIL（`MsdialExeNotFoundError` が `code` を受け取らない・`get_exe_path` が設定ファイルを読まない、など）

- [ ] **Step 3: `lipidmix/console/runner.py` を直す**

モジュール docstring の 3〜4 行目

```
exe_path は引数で注入できる（テスト用 fake の差し込みに使う）。
既定は環境変数 MSDIAL_EXE から取得し、未設定なら EnvironmentError を上げる。
```

を次に置き換える:

```
exe_path は引数で注入できる（テスト用 fake の差し込みに使う）。
既定は `lipidmix.core.user_config` の `msdial.exe`（環境変数 MSDIAL_EXE →
`lipidmix.local.toml` の `[msdial] exe`）から取得し、決められなければ
`MsdialExeNotFoundError`（EnvironmentError）を上げる。
```

import 群（`from pathlib import Path` の後）に追加:

```python
from lipidmix.core import user_config
```

`MsdialExeNotFoundError` と `get_exe_path` を次に置き換える:

```python
_EXE_WHAT = "MS-DIAL Console の実行体（MSDIALCUI.exe / MsdialConsoleApp.exe）"


class MsdialExeNotFoundError(EnvironmentError):
    """Console の実行体を決められない。

    `code` は封筒にそのまま載せる。未設定は `MSDIAL_EXE_NOT_FOUND`、設定ファイルが
    読めなければ `CONFIG_INVALID`。後者も同じ型で送るのは、呼び出し側の既存の
    `except EnvironmentError` をすべてそのまま生かすため。
    """

    def __init__(self, message: str, *, code: str = "MSDIAL_EXE_NOT_FOUND",
                 details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.details = details if details is not None else {}


def get_exe_path() -> str:
    """Console の実行体のパス（環境変数 MSDIAL_EXE → 設定ファイルの `[msdial] exe`）。"""
    try:
        setting = user_config.get_setting("msdial.exe")
    except user_config.ConfigInvalidError as exc:
        raise MsdialExeNotFoundError(exc.message, code=exc.code, details=exc.details()) from exc
    if setting is None:
        raise MsdialExeNotFoundError(
            "MS-DIAL Console の実行体が設定されていません。"
            + user_config.missing_hint("msdial.exe", _EXE_WHAT),
            details=user_config.describe_missing("msdial.exe"))
    return setting.value
```

- [ ] **Step 4: `lipidmix/tools/console_tools.py` を直す**

(a) `_msdial_exe_setup_help` を次に置き換える:

```python
def _msdial_exe_setup_help() -> dict:
    """Console の実行体が未設定の封筒に、誰が何をすべきかを機械可読で載せる。

    設定ファイル（`lipidmix.local.toml`）の編集も環境変数の設定も MCP クライアント
    からはできない。LLM に「設定してください」とだけ返すと、設定を試みて失敗するか
    黙って諦める。人間の作業であることを型で示し、コピペできる手順を渡す。
    設定ファイルは呼ばれるたびに読まれるので再起動は要らないが、環境変数は
    起動中のサーバに反映されない。
    """
    from lipidmix.console.runner import msdial_exe_candidates
    from lipidmix.core import user_config

    try:
        candidates = msdial_exe_candidates()
    except OSError:
        candidates = []
    config = user_config.config_file_path()
    if config.is_file():
        config_step = (f"{config} の [msdial] に exe = '<MSDIALCUI.exe のパス>' と書く"
                       "（単一引用符で囲む。サーバの再起動は不要）")
    else:
        config_step = (f"{user_config.EXAMPLE_FILENAME} を {config} として複製し、"
                       "[msdial] exe = '<MSDIALCUI.exe のパス>' と書く"
                       "（単一引用符で囲む。サーバの再起動は不要）")
    return {
        "human_action_required": True,
        "why": "設定ファイルの編集も環境変数の設定も MCP クライアントからはできません。",
        "candidates": candidates[:10],
        "how_to_set": [
            config_step,
            "PowerShell（恒久設定）: [Environment]::SetEnvironmentVariable("
            "'MSDIAL_EXE','<MSDIALCUI.exe のパス>','User')（設定後に MCP サーバの再起動が必要）",
        ],
        "restart_required": False,
        "restart_note": "設定ファイルなら次の呼び出しから効きます。環境変数で設定した場合だけ、"
                        "起動中のプロセスが変更を読み直さないので MCP サーバの再起動が必要です。",
    }
```

(b) `_msdial_exe_setup_help` の直前に 2 つのヘルパを足す:

```python
def _exe_error(exc: EnvironmentError) -> str:
    """`get_exe_path` の失敗を封筒にする。

    未設定（`MSDIAL_EXE_NOT_FOUND`）には人間がすべき手順を、設定ファイルが読めない
    （`CONFIG_INVALID`）ときはその行・列を載せる。
    """
    code = getattr(exc, "code", "MSDIAL_EXE_NOT_FOUND")
    details = dict(getattr(exc, "details", {}) or {})
    if code == "MSDIAL_EXE_NOT_FOUND":
        details = {**_msdial_exe_setup_help(), **details}
    return console_error(code, str(exc), details or None)


def _configured_exe() -> str:
    """封筒に載せるための設定済み実行体（決められなければ空文字）。"""
    from lipidmix.console import runner as console_runner
    try:
        return console_runner.get_exe_path()
    except EnvironmentError:
        return ""
```

(c) `console_plan` の

```python
    except EnvironmentError as exc:
        return console_error("MSDIAL_EXE_NOT_FOUND", str(exc),
                             _msdial_exe_setup_help())
```

を次に置き換える:

```python
    except EnvironmentError as exc:
        return _exe_error(exc)
```

(d) `console_run` の

```python
    except EnvironmentError as exc:
        update_status(resolved, "failed", error=str(exc))
        return console_error("MSDIAL_EXE_NOT_FOUND", str(exc))
    try:
        is_console = console_runner.is_console_exe(exe, raise_on_os_error=True)
```

を次に置き換える（2 つ目の `except OSError` 分岐はそのまま）:

```python
    except EnvironmentError as exc:
        update_status(resolved, "failed", error=str(exc))
        return _exe_error(exc)
    try:
        is_console = console_runner.is_console_exe(exe, raise_on_os_error=True)
```

(e) `console_method_template` の

```python
    exe = os.environ.get("MSDIAL_EXE") or None
```

を次に置き換える（直後の `resolve_lbm(..., env=os.environ)` は Task 3 で直すので触らない）:

```python
    from lipidmix.console import runner as console_runner
    try:
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        # exe はビルド生成物・exe フォルダからの LBM 推定にだけ使う。未設定でも
        # メソッドの宣言や [msdial] lbm で解決できれば足りるので、ここでは止めない。
        if getattr(exc, "code", "") == "CONFIG_INVALID":
            return _exe_error(exc)
        exe = None
```

(f) `_console_run_result` の `launch_failed` 分岐の

```python
            {**details, "exe": os.environ.get("MSDIAL_EXE", "")})
```

を次に置き換える:

```python
            {**details, "exe": _configured_exe()})
```

- [ ] **Step 5: `lipidmix/pipeline/service.py` の `_resolve_exe_path` を直す**

```python
def _resolve_exe_path() -> Path:
    from lipidmix.console import runner as console_runner
    try:
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        # 未設定は MSDIAL_EXE_NOT_FOUND、設定ファイルが読めなければ CONFIG_INVALID。
        # どちらも利用者が直す手掛かり（設定ファイルの場所・キー・行）を details に持つ。
        raise DomainError(getattr(exc, "code", "MSDIAL_EXE_NOT_FOUND"), str(exc),
                          dict(getattr(exc, "details", {}) or {})) from exc
    return Path(exe)
```

- [ ] **Step 6: 通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_console_config.py tests/test_msdial_exe_guidance.py tests/test_console_runner.py tests/test_console_method_template.py tests/test_pipeline_service.py -q`
Expected: 全て PASS。`test_msdial_exe_guidance.py` の 2 件（`human_action_required` と「再起動」の語、`how_to_set` が空でない）は新しい封筒でも満たされる。

- [ ] **Step 7: `console_tools.py` から不要になった `os` 参照を確認する**

Run: `grep -n "os\." lipidmix/tools/console_tools.py`
Expected: 残るのは `resolve_lbm(..., env=os.environ)` の 2 か所だけ（Task 3 で消える）。`import os` は Task 3 の最後に要否を判断するので、ここでは消さない。

- [ ] **Step 8: コミット（バックグラウンド）**

```bash
git add lipidmix/console/runner.py lipidmix/tools/console_tools.py lipidmix/pipeline/service.py tests/test_console_config.py
git commit -m "feat(console): Console の実行体を設定ファイルの [msdial] exe からも引く" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$SCRATCH/commit-t2.log" 2>&1; echo exit=$?; grep -E "passed|failed|FAILED" "$SCRATCH/commit-t2.log" | tail -3
```

---

### Task 3: 脂質ライブラリ（`[msdial] lbm`）

**Files:**
- Modify: `lipidmix/console/method_file.py`（import、`LbmResolution.source` のコメント、`resolve_lbm` の引数・docstring・`MSDIAL_LBM` 分岐・exe フォルダで見つからないときの文面）
- Modify: `lipidmix/tools/console_tools.py`（`console_plan` と `console_method_template` の `resolve_lbm` 呼び出し、新規ヘルパ `_lbm_setting`）
- Modify: `lipidmix/pipeline/inputs.py`（`_resolve_lbm_pinned`、import）
- Modify: `tests/test_console_method_file.py`（`resolve_lbm` の呼び出し全件）
- Modify: `tests/test_console_config.py`（テスト追加）

**Interfaces:**
- Consumes: Task 1 の `user_config.Setting` / `get_setting` / `setting_label` / `ConfigInvalidError`、Task 2 の `console_tools._exe_error`
- Produces:
  - `resolve_lbm(method_keys: dict[str, str], method_file: Path, omics: str, exe_path: str | None, lbm_setting: Setting | None = None, override: str | None = None) -> LbmResolution`（引数 `env` は廃止）
  - `LbmResolution.source` の値に `"config_file"` が加わる（`argument | method_file | build_tree | env | config_file | exe_dir | not_required`）
  - `console_tools._lbm_setting() -> Setting | None | str`（`str` は封筒）

- [ ] **Step 1: 既存テストの呼び出しを新しい引数に直す**

`tests/test_console_method_file.py` で:

1. `env={}` を渡しているすべての呼び出しから、その引数を消す。改行位置に注意して、次の 2 通りを置換する:
   - `, env={})` → `)`
   - `exe_path=None, env={})` → `exe_path=None)`（上の置換で済むが、残っていないことを grep で確認）
2. `env={"MSDIAL_LBM": X}` を渡している 5 件（`test_resolve_lbm_falls_back_to_env`、`test_resolve_lbm_errors_when_env_path_is_missing`、`test_resolve_lbm_not_required_for_metabolomics`、`test_resolve_lbm_prefers_the_build_tree_over_the_installed_env_path`、`test_resolve_lbm_build_tree_excludes_msdial4`）は `lbm_setting=_env_lbm(X)` にする。
3. ファイル冒頭の import 群に追加:

```python
from lipidmix.core.user_config import Setting


def _env_lbm(value) -> Setting:
    """環境変数 MSDIAL_LBM から来た値（resolve_lbm は出どころごと受け取る）。"""
    return Setting(key="msdial.lbm", value=str(value), source="env",
                   env_var="MSDIAL_LBM", config_file=None)
```

確認: `grep -n "env=" tests/test_console_method_file.py` が何も返さないこと。

- [ ] **Step 2: 失敗するテストを足す**

`tests/test_console_method_file.py` の `test_resolve_lbm_errors_when_env_path_is_missing` の直後に:

```python
def test_resolve_lbm_takes_the_value_from_the_config_file(tmp_path):
    lib = tmp_path / "conf.lbm2"
    lib.touch()
    exe = _exe_dir_with(tmp_path)
    setting = Setting(key="msdial.lbm", value=str(lib), source="config_file",
                      env_var="MSDIAL_LBM", config_file=str(tmp_path / "lipidmix.local.toml"))
    res = resolve_lbm({"lbm file path": ""}, tmp_path / "param.txt",
                      omics="lipidomics", exe_path=str(exe), lbm_setting=setting)
    assert res.source == "config_file"
    assert Path(res.path) == lib


def test_resolve_lbm_names_the_config_key_when_its_file_is_missing(tmp_path):
    exe = _exe_dir_with(tmp_path, "one.lbm2")
    setting = Setting(key="msdial.lbm", value=str(tmp_path / "gone.lbm2"), source="config_file",
                      env_var="MSDIAL_LBM", config_file=str(tmp_path / "lipidmix.local.toml"))
    res = resolve_lbm({}, tmp_path / "param.txt", omics="lipidomics",
                      exe_path=str(exe), lbm_setting=setting)
    assert res.error_code == "LBM_NOT_FOUND"
    assert "lipidmix.local.toml の [msdial] lbm" in res.message


def test_resolve_lbm_names_the_env_var_when_its_file_is_missing(tmp_path):
    exe = _exe_dir_with(tmp_path, "one.lbm2")
    res = resolve_lbm({}, tmp_path / "param.txt", omics="lipidomics",
                      exe_path=str(exe), lbm_setting=_env_lbm(tmp_path / "gone.lbm2"))
    assert "環境変数 MSDIAL_LBM" in res.message
```

`tests/test_console_config.py` の末尾に:

```python
# ---------- LBM（[msdial] lbm） ----------

def test_template_uses_lbm_from_the_config_file_while_exe_comes_from_env(tmp_path, cfg, monkeypatch):
    """出どころが混ざる構成: exe は環境変数、LBM は設定ファイル。"""
    from lipidmix.tools.console_tools import console_method_template
    app = tmp_path / "app"
    app.mkdir()
    (app / "bundled.lbm2").touch()  # exe フォルダの LBM より設定ファイルの LBM が先
    lib = tmp_path / "libs" / "chosen.lbm2"
    lib.parent.mkdir()
    lib.touch()
    monkeypatch.setenv("MSDIAL_EXE", str(app / "MSDIALCUI.exe"))
    cfg.write_text(f"[msdial]\nlbm = '{lib}'\n", encoding="utf-8")
    out = tmp_path / "param_POS.txt"
    parsed = json.loads(console_method_template(
        out_path=str(out), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["lbm"]["source"] == "config_file"
    assert "chosen.lbm2" in out.read_text(encoding="ascii")


def test_template_reports_config_invalid_from_the_lbm_lookup(tmp_path, cfg, monkeypatch):
    """exe を環境変数で渡していても、LBM を引くときに壊れた設定ファイルに当たれば止める。"""
    from lipidmix.tools.console_tools import console_method_template
    app = tmp_path / "app"
    app.mkdir()
    monkeypatch.setenv("MSDIAL_EXE", str(app / "MSDIALCUI.exe"))
    cfg.write_text("[msdial\n", encoding="utf-8")
    parsed = json.loads(console_method_template(
        out_path=str(tmp_path / "param_POS.txt"), polarity="positive",
        based_on=str(_neg_param(tmp_path / "neg_param_1.txt"))))
    assert parsed["error"]["code"] == "CONFIG_INVALID"


def test_pipeline_pins_the_lbm_from_the_config_file(tmp_path, cfg):
    from lipidmix.pipeline import inputs
    lib = tmp_path / "chosen.lbm2"
    lib.write_bytes(b"lbm")
    cfg.write_text(f"[msdial]\nlbm = '{lib}'\n", encoding="utf-8")
    pinned = inputs._resolve_lbm_pinned({}, tmp_path / "param.txt",
                                        str(tmp_path / "app" / "MSDIALCUI.exe"), {}, tmp_path)
    assert pinned["source"] == "config_file"
    assert Path(pinned["path"]) == lib.resolve()


def test_pipeline_lbm_lookup_wraps_config_invalid(tmp_path, cfg):
    from lipidmix.core.atomic_io import DomainError
    from lipidmix.pipeline import inputs
    cfg.write_text("[msdial\n", encoding="utf-8")
    with pytest.raises(DomainError) as exc:
        inputs._resolve_lbm_pinned({}, tmp_path / "param.txt",
                                   str(tmp_path / "app" / "MSDIALCUI.exe"), {}, tmp_path)
    assert exc.value.code == "CONFIG_INVALID"
```

- [ ] **Step 3: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_console_method_file.py tests/test_console_config.py -q`
Expected: FAIL（`resolve_lbm() got an unexpected keyword argument 'lbm_setting'` など）

- [ ] **Step 4: `lipidmix/console/method_file.py` を直す**

import 群（`from lipidmix.core.atomic_io import DomainError` の次）に追加:

```python
from lipidmix.core.user_config import Setting, setting_label
```

`LbmResolution` の

```python
    source: str  # argument | method_file | build_tree | env | exe_dir | not_required
```

を

```python
    source: str  # argument | method_file | build_tree | env | config_file | exe_dir | not_required
```

にする。

`resolve_lbm` のシグネチャと docstring 冒頭を次に置き換える（docstring の 2 段落目以降はそのまま）:

```python
def resolve_lbm(
    method_keys: dict[str, str],
    method_file: Path,
    omics: str,
    exe_path: str | None,
    lbm_setting: Setting | None = None,
    override: str | None = None,
) -> LbmResolution:
    """脂質ライブラリのパスを GUI と同じ規則で解決する。

    順に: 明示引数 → メソッドファイルの宣言 → **ビルド生成物** → `lbm_setting`
    （環境変数 MSDIAL_LBM → 設定ファイルの `[msdial] lbm`。呼び出し側が
    `user_config.get_setting("msdial.lbm")` で引いて渡す）→ Console の実行体と同じ
    フォルダ。ビルド生成物を `lbm_setting` より上に置くのは、
    この環境の Console がソースからのビルドで、ライブラリもそのツリー内の
    新しいものを使うため（インストール版より優先する）。
```

`MSDIAL_LBM` 分岐

```python
    from_env = (env.get("MSDIAL_LBM") or "").strip()
    if from_env:
        candidate = Path(from_env)
        if candidate.is_file():
            return LbmResolution(path=_absolute(candidate), source="env")
        return LbmResolution(
            path=None, source="env", error_code="LBM_NOT_FOUND",
            message=f"環境変数 MSDIAL_LBM が指すファイルがありません: {from_env}")
```

を次に置き換える:

```python
    if lbm_setting is not None:
        candidate = Path(lbm_setting.value)
        if candidate.is_file():
            return LbmResolution(path=_absolute(candidate), source=lbm_setting.source)
        return LbmResolution(
            path=None, source=lbm_setting.source, error_code="LBM_NOT_FOUND",
            message=f"{setting_label(lbm_setting)} が指すファイルがありません: {lbm_setting.value}")
```

exe フォルダに見つからないときの文面の末尾

```python
                "MS-DIAL のインストールフォルダにある .lbm2 のパスを環境変数"
                "MSDIAL_LBM に設定するか、lbm_file 引数で渡してください。"))
```

を

```python
                "MS-DIAL のインストールフォルダにある .lbm2 のパスを設定ファイル"
                "（lipidmix.local.toml）の [msdial] lbm か環境変数 MSDIAL_LBM に"
                "設定するか、lbm_file 引数で渡してください。"))
```

にする。exe 未設定の文面 `"脂質ライブラリ（.lbm2）を解決できません。MSDIAL_EXE が未設定です。"` は `"脂質ライブラリ（.lbm2）を解決できません。Console の実行体（[msdial] exe / MSDIAL_EXE）が未設定です。"` にする。

- [ ] **Step 5: `lipidmix/tools/console_tools.py` の呼び出しを直す**

`_exe_error` の直後にヘルパを足す:

```python
def _lbm_setting():
    """`[msdial] lbm` / MSDIAL_LBM の設定。設定ファイルが読めなければ封筒（str）を返す。"""
    from lipidmix.core import user_config
    try:
        return user_config.get_setting("msdial.lbm")
    except user_config.ConfigInvalidError as exc:
        return console_error(exc.code, exc.message, exc.details())
```

`console_plan` の

```python
    lbm = method_file_mod.resolve_lbm(
        method_keys, mf, omics=omics, exe_path=exe, env=os.environ, override=lbm_file)
```

を

```python
    lbm_setting = _lbm_setting()
    if isinstance(lbm_setting, str):
        return lbm_setting
    lbm = method_file_mod.resolve_lbm(
        method_keys, mf, omics=omics, exe_path=exe, lbm_setting=lbm_setting, override=lbm_file)
```

に、`console_method_template` の

```python
    lbm = method_file_mod.resolve_lbm(
        src_keys, src, omics=omics, exe_path=exe, env=os.environ)
```

を

```python
    lbm_setting = _lbm_setting()
    if isinstance(lbm_setting, str):
        return lbm_setting
    lbm = method_file_mod.resolve_lbm(
        src_keys, src, omics=omics, exe_path=exe, lbm_setting=lbm_setting)
```

にする。その後 `grep -n "os\." lipidmix/tools/console_tools.py` が何も返さなければ、先頭の `import os` を消す（他に `os` を使っていれば残す）。

- [ ] **Step 6: `lipidmix/pipeline/inputs.py` を直す**

import 群の `from lipidmix.core import app_control` の次に追加:

```python
from lipidmix.core import user_config
```

`_resolve_lbm_pinned` の

```python
    lbm = method_file_mod.resolve_lbm(
        method_keys, method_path, omics=_OMICS, exe_path=exe_path,
        env=dict(os.environ), override=override)
```

を

```python
    try:
        lbm_setting = user_config.get_setting("msdial.lbm")
    except user_config.ConfigInvalidError as exc:
        raise DomainError(exc.code, exc.message, exc.details()) from exc
    lbm = method_file_mod.resolve_lbm(
        method_keys, method_path, omics=_OMICS, exe_path=exe_path,
        lbm_setting=lbm_setting, override=override)
```

にする。`import os` は他で使っているので残す。

- [ ] **Step 7: 通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_console_method_file.py tests/test_console_config.py tests/test_console_plan_method.py tests/test_console_method_template.py tests/test_pipeline_inputs.py tests/test_lcms_profile_inputs.py tests/test_console_runner.py -q`
Expected: 全て PASS

Run: `grep -rn "resolve_lbm(" lipidmix tests | grep "env="`
Expected: 何も出ない

- [ ] **Step 8: コミット（バックグラウンド）**

```bash
git add lipidmix/console/method_file.py lipidmix/tools/console_tools.py lipidmix/pipeline/inputs.py tests/test_console_method_file.py tests/test_console_config.py
git commit -m "feat(console): 脂質ライブラリを設定ファイルの [msdial] lbm からも引く" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$SCRATCH/commit-t3.log" 2>&1; echo exit=$?; grep -E "passed|failed|FAILED" "$SCRATCH/commit-t3.log" | tail -3
```

---

### Task 4: 研究室の参照ライブラリ（`[library] msp_positive` / `msp_negative`）

**Files:**
- Modify: `lipidmix/core/path_resolvers.py`（`LIBRARY_ENV_VARS` 周辺、`LibraryPathError`、`_library_from_env` → `_library_from_setting`、`resolve_library_path` の docstring と両極性判定）
- Modify: `lipidmix/library/tools.py`（`library_load` の docstring とエラー戻り値）
- Modify: `lipidmix/library/store.py`（CLI の `--ion-mode` の help）
- Create: `tests/test_library_config.py`

**Interfaces:**
- Consumes: Task 1 の `user_config.get_setting` / `setting_label` / `describe_missing` / `missing_hint` / `SETTINGS` / `ConfigInvalidError`
- Produces:
  - `path_resolvers.LIBRARY_SETTING_KEYS = {"positive": "library.msp_positive", "negative": "library.msp_negative"}`
  - `path_resolvers.LIBRARY_ENV_VARS`（同名で残す。値は `user_config.SETTINGS` から導く）
  - `LibraryPathError(code: str, message: str, details: dict | None = None)`、属性 `details: dict`
  - `library_load` のエラー戻り値に `details`（あれば）が加わる: `{"status": "error", "code": ..., "message": ..., "details": {...}}`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_library_config.py`:

```python
"""研究室の参照ライブラリを設定ファイルから引く（spec 2026-10-06 §6・§7）。

置き場所（ディレクトリ）は戻り値に出さない規則を、設定ファイル経由でも守る。
"""
from __future__ import annotations

import json

import pytest

from lipidmix.core.path_resolvers import LibraryPathError, resolve_library_path

BS = "\\"


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    path = tmp_path / "conf" / "lipidmix.local.toml"
    path.parent.mkdir()
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(path))
    # data ディレクトリの *.msp / *.dbs を拾わせない
    from lipidmix.core import mcp_core
    empty = tmp_path / "data"
    empty.mkdir()
    monkeypatch.setattr(mcp_core, "DATA_DIR", empty)
    return path


def _msp(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("NAME: x\nPRECURSORMZ: 100\nNum Peaks: 0\n\n", encoding="utf-8")
    return path


def test_ion_mode_reads_the_config_file(cfg, tmp_path):
    pos = _msp(tmp_path / "lab" / "pos.msp")
    cfg.write_text(f"[library]\nmsp_positive = '{pos}'\n", encoding="utf-8")
    assert resolve_library_path(ion_mode="positive") == str(pos)


def test_single_configured_polarity_is_used_without_ion_mode(cfg, tmp_path):
    neg = _msp(tmp_path / "lab" / "neg.msp")
    cfg.write_text(f"[library]\nmsp_negative = '{neg}'\n", encoding="utf-8")
    assert resolve_library_path() == str(neg)


def test_both_polarities_across_sources_are_ambiguous(cfg, tmp_path, monkeypatch):
    """片方が環境変数・片方が設定ファイルでも、両方あれば ion_mode を求める。"""
    pos = _msp(tmp_path / "lab" / "pos.msp")
    neg = _msp(tmp_path / "lab" / "neg.msp")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(pos))
    cfg.write_text(f"[library]\nmsp_negative = '{neg}'\n", encoding="utf-8")
    with pytest.raises(LibraryPathError) as exc:
        resolve_library_path()
    assert exc.value.code == "MSP_AMBIGUOUS"
    assert "環境変数 MSDIAL_MSP_POS" in exc.value.message
    assert "[library] msp_negative" in exc.value.message


def test_missing_file_from_the_config_names_the_key_not_the_directory(cfg, tmp_path):
    gone = tmp_path / "secret-lab-share" / "neg.msp"
    cfg.write_text(f"[library]\nmsp_negative = '{gone}'\n", encoding="utf-8")
    with pytest.raises(LibraryPathError) as exc:
        resolve_library_path(ion_mode="negative")
    assert exc.value.code == "MSP_ENV_NOT_FOUND"
    assert "lipidmix.local.toml の [library] msp_negative" in exc.value.message
    assert "neg.msp" in exc.value.message
    dumped = exc.value.message + json.dumps(exc.value.details, ensure_ascii=False)
    assert "secret-lab-share" not in dumped


def test_broken_config_file_is_config_invalid(cfg):
    cfg.write_text("[library\n", encoding="utf-8")
    with pytest.raises(LibraryPathError) as exc:
        resolve_library_path(ion_mode="positive")
    assert exc.value.code == "CONFIG_INVALID"


def test_library_load_missing_config_msp_does_not_leak_the_directory(cfg, tmp_path):
    from lipidmix.library.tools import library_load
    gone = tmp_path / "secret-lab-share" / "pos.msp"
    cfg.write_text(f"[library]\nmsp_positive = '{gone}'\n", encoding="utf-8")
    raw = library_load(ion_mode="positive")
    payload = json.loads(raw)
    assert payload["code"] == "MSP_ENV_NOT_FOUND"
    assert payload["details"]["setting"] == "library.msp_positive"
    assert "secret-lab-share" not in raw


def test_library_load_with_nothing_configured_says_where_to_write(cfg):
    from lipidmix.library.tools import library_load
    payload = json.loads(library_load(ion_mode="negative"))
    assert payload["status"] == "error"
    assert "[library] msp_negative" in payload["message"]
    assert payload["details"]["setting"] == "library.msp_negative"
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_library_config.py -q`
Expected: FAIL（設定ファイルを読まない・`details` が無い）

注意: `cfg` fixture は `mcp_core.DATA_DIR` を差し替えている。`list_data_files` は `mcp_core.DATA_DIR` を module 修飾で動的参照している（path_resolvers.py の `target_dir = ... else mcp_core.DATA_DIR`）ので、これで data ディレクトリの探索は空になる。

- [ ] **Step 3: `lipidmix/core/path_resolvers.py` を直す**

import 群の `from lipidmix.core import mcp_core` の次に追加:

```python
from lipidmix.core import user_config
```

モジュール docstring 3 行目の `依存は mcp_core（DATA_DIR）と stdlib のみの下位レイヤ。` を `依存は mcp_core（DATA_DIR）・user_config（外部資産の場所）と stdlib のみの下位レイヤ。` にする。

`LIBRARY_ENV_VARS` の定義とその上のコメント、`LibraryPathError`、`_library_from_env` を次に置き換える:

```python
#: 極性 → 研究室参照ライブラリ（`.msp`）の設定キー。値は環境変数か
#: `lipidmix.local.toml` の `[library]` で指す（`lipidmix.core.user_config`）。
#: ライブラリ本体はリポジトリの外（外部流出禁止の資産）に置き、場所だけを教える。
LIBRARY_SETTING_KEYS = {"positive": "library.msp_positive", "negative": "library.msp_negative"}
#: 同じ値を指す環境変数（正準は user_config.SETTINGS。既存の参照のため名前を残す）。
LIBRARY_ENV_VARS = {mode: user_config.SETTINGS[key] for mode, key in LIBRARY_SETTING_KEYS.items()}

_MSP_WHAT = "研究室の参照ライブラリ（.msp）"


class LibraryPathError(Exception):
    """参照ライブラリを 1 つに決められない（または指定先が無い）。

    `code` は機械可読（`LIBRARY_NOT_FOUND` / `MSP_ENV_NOT_FOUND` /
    `MSP_AMBIGUOUS` / `INVALID_ION_MODE` / `CONFIG_INVALID`）。`message` と `details`
    にはファイル名・設定キー・環境変数名・設定ファイルのパスだけを載せ、ライブラリの
    置き場所（ディレクトリ）は載せない——戻り値は LLM の文脈に入る。
    `MSP_ENV_NOT_FOUND` は値が設定ファイルから来た場合も同じコード（互換のため据え置き）。
    """

    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details if details is not None else {}


def _library_setting(ion_mode: str) -> user_config.Setting | None:
    """極性の設定（環境変数 → 設定ファイル）。設定ファイルが読めなければ CONFIG_INVALID。"""
    try:
        return user_config.get_setting(LIBRARY_SETTING_KEYS[ion_mode])
    except user_config.ConfigInvalidError as exc:
        raise LibraryPathError(exc.code, exc.message, exc.details()) from exc


def _library_from_setting(ion_mode: str) -> str | None:
    """極性の設定が指すパス。未設定なら None、指す先が無ければ例外。"""
    setting = _library_setting(ion_mode)
    if setting is None:
        return None
    if not os.path.isfile(setting.value):
        label = user_config.setting_label(setting)
        raise LibraryPathError(
            "MSP_ENV_NOT_FOUND",
            f"{label} が指すファイル（{os.path.basename(setting.value)}）がありません。"
            f"{label} を確認してください。",
            user_config.describe_missing(LIBRARY_SETTING_KEYS[ion_mode], setting),
        )
    return setting.value
```

`resolve_library_path` の docstring の 2・4 番

```
    2. `ion_mode` を指定したら、その極性の環境変数（`LIBRARY_ENV_VARS`）。
       未設定なら 3 以降へ落ちる。
    3. data ディレクトリの `*_Loaded.msp2.dbs`（その run が実際に使った参照）。
    4. 極性の環境変数。1 つだけ設定されていればそれ、両方なら `MSP_AMBIGUOUS`
       （`ion_mode` の指定を求める）。
```

を次に置き換える:

```
    2. `ion_mode` を指定したら、その極性の設定（環境変数 `MSDIAL_MSP_POS` /
       `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の `[library] msp_positive` /
       `msp_negative`。`LIBRARY_SETTING_KEYS`）。未設定なら 3 以降へ落ちる。
    3. data ディレクトリの `*_Loaded.msp2.dbs`（その run が実際に使った参照）。
    4. 極性の設定。1 つだけ設定されていればそれ、両方なら `MSP_AMBIGUOUS`
       （`ion_mode` の指定を求める）。出どころ（環境変数か設定ファイルか）は問わない。
```

本体の

```python
    if ion_mode is not None:
        from_env = _library_from_env(ion_mode)
        if from_env:
            return from_env
```

を

```python
    if ion_mode is not None:
        from_setting = _library_from_setting(ion_mode)
        if from_setting:
            return from_setting
```

に、

```python
    configured = [mode for mode, var in LIBRARY_ENV_VARS.items() if (os.environ.get(var) or "").strip()]
    if len(configured) > 1:
        names = " / ".join(LIBRARY_ENV_VARS[mode] for mode in configured)
        raise LibraryPathError(
            "MSP_AMBIGUOUS",
            f"参照ライブラリが極性ごとに設定されています（{names}）。"
            f"ion_mode（{' / '.join(configured)}）を指定してください。",
        )
    if configured:
        return _library_from_env(configured[0])
```

を

```python
    settings = {mode: _library_setting(mode) for mode in LIBRARY_SETTING_KEYS}
    configured = [mode for mode, setting in settings.items() if setting is not None]
    if len(configured) > 1:
        names = " / ".join(user_config.setting_label(settings[mode]) for mode in configured)
        raise LibraryPathError(
            "MSP_AMBIGUOUS",
            f"参照ライブラリが極性ごとに設定されています（{names}）。"
            f"ion_mode（{' / '.join(configured)}）を指定してください。",
        )
    if configured:
        return _library_from_setting(configured[0])
```

に置き換える。

- [ ] **Step 4: `lipidmix/library/tools.py` の `library_load` を直す**

docstring の解決順の段落

```
    解決順: `file_path` の明示 → `ion_mode`（`"positive"` / `"negative"`）に対応する
    環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → データディレクトリの
    `*_Loaded.msp2.dbs` → 設定済みの環境変数（両方あれば `ion_mode` を求める）→
    データディレクトリの `*.msp`。
```

を

```
    解決順: `file_path` の明示 → `ion_mode`（`"positive"` / `"negative"`）に対応する
    設定（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の
    `[library] msp_positive` / `msp_negative`）→ データディレクトリの
    `*_Loaded.msp2.dbs` → 設定済みの極性（両方あれば `ion_mode` を求める）→
    データディレクトリの `*.msp`。
```

にする。

エラー処理

```python
    except LibraryPathError as exc:
        return json_payload({"status": "error", "code": exc.code, "message": exc.message})
    if not resolved:
        return json_payload({
            "status": "error",
            "message": "データディレクトリに参照ライブラリ（*_Loaded.msp2.dbs または *.msp）が見つかりませんでした。",
        })
```

を次に置き換える:

```python
    except LibraryPathError as exc:
        payload = {"status": "error", "code": exc.code, "message": exc.message}
        if exc.details:
            payload["details"] = exc.details
        return json_payload(payload)
    if not resolved:
        key = LIBRARY_SETTING_KEYS.get(ion_mode or "positive", "library.msp_positive")
        hint = user_config.missing_hint(key, "研究室の参照ライブラリ（.msp）")
        if ion_mode is None:
            hint += " 負イオンは [library] msp_negative です。"
        return json_payload({
            "status": "error",
            "message": "参照ライブラリが見つかりませんでした（データディレクトリにも "
                       "*_Loaded.msp2.dbs / *.msp がありません）。file_path で指定するか、" + hint,
            "details": user_config.describe_missing(key),
        })
```

import 行

```python
from lipidmix.core.path_resolvers import LibraryPathError, resolve_dcl_file_path, resolve_library_path
```

を

```python
from lipidmix.core import user_config
from lipidmix.core.path_resolvers import (
    LIBRARY_SETTING_KEYS,
    LibraryPathError,
    resolve_dcl_file_path,
    resolve_library_path,
)
```

にする。`ion_mode` が `"POSITIVE"` のような大文字で来ても `resolve_library_path` が正規化済みなので、ここは `.get(..., 既定)` で安全に落ちる。

- [ ] **Step 5: `lipidmix/library/store.py` の CLI help を直す**

```python
                        help="環境変数 MSDIAL_MSP_POS / MSDIAL_MSP_NEG のどちらを使うか")
    parser.add_argument("--file", help="ライブラリのパス（環境変数より優先）")
```

を

```python
                        help="どちらの極性の設定（MSDIAL_MSP_POS / MSDIAL_MSP_NEG か "
                             "lipidmix.local.toml の [library] msp_positive / msp_negative）を使うか")
    parser.add_argument("--file", help="ライブラリのパス（極性の設定より優先）")
```

にする。

- [ ] **Step 6: 通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_library_config.py tests/test_library_session.py tests/test_library_tools.py tests/test_library_store.py -q`
Expected: 全て PASS。既存の `test_library_session.py:103`（`"MSDIAL_MSP_NEG" in exc.value.message`）は「環境変数 MSDIAL_MSP_NEG が指す…」で満たされる。

- [ ] **Step 7: コミット（バックグラウンド）**

```bash
git add lipidmix/core/path_resolvers.py lipidmix/library/tools.py lipidmix/library/store.py tests/test_library_config.py
git commit -m "feat(library): 研究室の参照ライブラリを設定ファイルの [library] からも引く" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$SCRATCH/commit-t4.log" 2>&1; echo exit=$?; grep -E "passed|failed|FAILED" "$SCRATCH/commit-t4.log" | tail -3
```

---

### Task 5: 雛形・worktree スクリプト・文書・記録

**Files:**
- Create: `lipidmix.example.toml`
- Delete: `.env.example`
- Create: `tests/test_example_config.py`
- Modify: `scripts/new-worktree.sh`
- Modify: `README.md`（Quick start）
- Modify: `USAGE.md`（`console_method_template` 行、`library_load` 行）
- Modify: `docs/output_format/library.md`（エラー戻り値の段落）
- Modify: `docs/workflow/library.md`（`library_load` の前提）
- Modify: `docs/cli.md`（`--ion-mode` の行）
- Modify: `CLAUDE.md`（環境変数の段落、研究室ライブラリの段落、stdlib leaf の一覧）
- 追跡外: main ツリーの `docs/HISTRY.md` / `docs/task.md`、vault の流れ図

**Interfaces:**
- Consumes: Task 1〜4 の全て
- Produces: なし（利用者向けの入口）

- [ ] **Step 1: 雛形が実際に読めることを縛るテストを書く**

`tests/test_example_config.py`:

```python
"""雛形 lipidmix.example.toml がそのまま設定ファイルとして読めること。

雛形が壊れていると、外部の利用者は複製した最初の一歩で CONFIG_INVALID に当たる。
"""
from __future__ import annotations

import shutil

from lipidmix.core import user_config


def test_example_is_a_valid_config_with_only_known_keys(tmp_path, monkeypatch):
    copied = tmp_path / "lipidmix.local.toml"
    shutil.copyfile(user_config.REPO_ROOT / user_config.EXAMPLE_FILENAME, copied)
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(copied))
    for key in user_config.SETTINGS:
        user_config.get_setting(key)  # CONFIG_INVALID を送出しない
    assert "unknown_keys" not in user_config.describe_missing("msdial.exe")
    assert user_config.get_setting("msdial.exe").value.endswith("MSDIALCUI.exe")


def test_env_example_is_gone():
    """誰にも読まれず LIPIDMIX_DATA_DIR を勧めていた .env.example は雛形に置き換えた。"""
    assert not (user_config.REPO_ROOT / ".env.example").exists()
```

Run: `C:/Python314/python.exe -m pytest tests/test_example_config.py -q`
Expected: FAIL（雛形が無い）

- [ ] **Step 2: 雛形を作り、`.env.example` を消す**

`lipidmix.example.toml`:

```toml
# ms-data-parser の外部資産の場所。
#
# このファイルを lipidmix.local.toml という名前でリポジトリ直下（server.py の隣）に
# 複製して編集する。lipidmix.local.toml は git に追跡されない。
#
# - Windows のパスは単一引用符で囲む（二重引用符だと \ がエスケープとして読まれる）。
# - 相対パスは、このファイルのあるフォルダ基準で解く。
# - 環境変数（MSDIAL_EXE / MSDIAL_LBM / MSDIAL_MSP_POS / MSDIAL_MSP_NEG）が
#   設定されていれば、そちらが優先される。
# - 書き換えは次のツール呼び出しから効く（MCP サーバの再起動は不要）。

[msdial]
# MS-DIAL Console の実行体。console_* / pipeline_* ツールを使うときだけ必要。
exe = 'C:\MS-DIAL\MSDIALCUI.exe'

# 脂質ライブラリ（.lbm2）。任意。書かなければメソッドファイルの宣言や、
# Console の実行体と同じフォルダの .lbm2 から推定する。
# lbm = 'C:\MS-DIAL\Lipids.lbm2'

[library]
# library_load(ion_mode=...) が読む参照ライブラリ（.msp）。極性ごとに 1 つ。
# msp_positive = 'D:\library\pos.msp'
# msp_negative = 'D:\library\neg.msp'
```

```bash
git rm .env.example
```

Run: `C:/Python314/python.exe -m pytest tests/test_example_config.py -q`
Expected: PASS

- [ ] **Step 3: `scripts/new-worktree.sh` で main の設定ファイルを worktree に届ける**

`main_w=$(win "$main")` と `wt_w=$(win "$path")` の 2 行の直後に追加:

```sh
# 外部資産の場所の設定（lipidmix.local.toml、追跡外）は worktree に来ない。
# main ツリーにあれば LIPIDMIX_CONFIG でそれを指す（設定を分裂させない）。
config_env=""
if [ -f "$main/lipidmix.local.toml" ]; then
    config_env="\"LIPIDMIX_CONFIG\": \"$main_w/lipidmix.local.toml\","
fi
```

heredoc 内の `"env": {` の次の行（`"LIPIDMIX_DATA_DIR": ...` の前）に 1 行足す:

```
        $config_env
```

（`config_env` が空なら空行になるだけで JSON は壊れない。）

確認:

```bash
sh -n scripts/new-worktree.sh && echo syntax-ok
```

Expected: `syntax-ok`。さらに main ツリーに `lipidmix.local.toml` が無い状態と、一時的に空ファイルを置いた状態の両方で、生成部分だけを抜き出して確かめる:

```bash
main=$(pwd); main_w=$(cygpath -m "$main"); config_env=""; [ -f "$main/lipidmix.local.toml" ] && config_env="\"LIPIDMIX_CONFIG\": \"$main_w/lipidmix.local.toml\","; printf '{"env": {\n        %s\n        "A": "b"}}\n' "$config_env" | C:/Python314/python.exe -c "import json,sys; print(json.load(sys.stdin))"
```

Expected: どちらの状態でも JSON として読める（置いた空ファイルは確認後に消す）。

- [ ] **Step 4: README の Quick start に設定ファイルの段落を足す**

`README.md` の `To run the parsers without an MCP client, see [docs/cli.md](docs/cli.md).` の直前に追加:

```markdown
### MS-DIAL Console and reference libraries (optional)

Running MS-DIAL Console from raw data (`console_*`, `pipeline_*`) or matching against your own
`.msp` library (`library_load`) needs to know where those files live. Copy
`lipidmix.example.toml` to `lipidmix.local.toml` next to `server.py` and fill in the paths
(Windows paths in single quotes). The file is not tracked by git and is read on every call, so
edits take effect without restarting the server. The environment variables `MSDIAL_EXE`,
`MSDIAL_LBM`, `MSDIAL_MSP_POS` and `MSDIAL_MSP_NEG` still work and take precedence over the file.
If something is missing, the tool's error names the file and the key to fill in.
```

- [ ] **Step 5: USAGE.md を直す**

`console_method_template` 行の `(解決順は `console_plan` と同じで、ビルド生成物が `MSDIAL_LBM` より優先)` を
`(解決順は `console_plan` と同じで、ビルド生成物が `[msdial] lbm` / `MSDIAL_LBM` より優先)` にする。

`library_load` 行の
`` `ion_mode`(`"positive"`/`"negative"`)を渡すと環境変数 `MSDIAL_MSP_POS`/`MSDIAL_MSP_NEG` が指す研究室ライブラリを読む。解決順は `file_path` → `ion_mode` の環境変数 → `*_Loaded.msp2.dbs` → 設定済みの環境変数 → `*.msp` で、 ``
を
`` `ion_mode`(`"positive"`/`"negative"`)を渡すと、その極性の設定(環境変数 `MSDIAL_MSP_POS`/`MSDIAL_MSP_NEG`、無ければ `lipidmix.local.toml` の `[library] msp_positive`/`msp_negative`)が指す研究室ライブラリを読む。解決順は `file_path` → `ion_mode` の設定 → `*_Loaded.msp2.dbs` → 設定済みの極性 → `*.msp` で、 ``
にする。

`console_plan` の行の末尾（行内の最後の `。` の後）に次の文を足す:
`` Console の実行体は環境変数 `MSDIAL_EXE` か `lipidmix.local.toml` の `[msdial] exe`(雛形 `lipidmix.example.toml`)で指す。未設定なら `MSDIAL_EXE_NOT_FOUND` の `details` に設定ファイルの場所・キー・手順が入り、設定ファイルが読めなければ `CONFIG_INVALID`(行・列つき)を返す。 ``

Run: `C:/Python314/python.exe -m pytest tests/test_readme_links.py -q`
Expected: PASS

- [ ] **Step 6: `docs/output_format/library.md` と `docs/workflow/library.md` と `docs/cli.md` を直す**

`docs/output_format/library.md` のエラー戻り値の段落

```
**エラー戻り値**: ライブラリを 1 つに決められないときは `{"status": "error", "code": ..., "message": ...}`
を返し、store は差し替えない。`code` は `MSP_AMBIGUOUS`（候補が複数。`ion_mode` か
`file_path` を指定する）/ `MSP_ENV_NOT_FOUND`（環境変数 `MSDIAL_MSP_POS` /
`MSDIAL_MSP_NEG` の指す先が無い）/ `LIBRARY_NOT_FOUND`（明示した `file_path` が無い）/
`INVALID_ION_MODE`。`message` にはファイル名と環境変数名だけが入り、置き場所は入らない。
```

を次に置き換える:

```
**エラー戻り値**: ライブラリを 1 つに決められないときは `{"status": "error", "code": ..., "message": ..., "details": ...}`
を返し、store は差し替えない。`code` は `MSP_AMBIGUOUS`（候補が複数。`ion_mode` か
`file_path` を指定する）/ `MSP_ENV_NOT_FOUND`（極性の設定——環境変数 `MSDIAL_MSP_POS` /
`MSDIAL_MSP_NEG` か `lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`——の
指す先が無い。設定ファイル由来でも同じコード）/ `LIBRARY_NOT_FOUND`（明示した `file_path` が無い）/
`INVALID_ION_MODE` / `CONFIG_INVALID`（設定ファイルが読めない。`details` に `config_file`
`line` `column`）。`details` は `MSP_ENV_NOT_FOUND` と、何も見つからなかったときに付き、
`setting`（設定キー）`env_var` `source`（`env` / `config_file` / null）`config_file`
`config_file_exists` `example`、あれば `unknown_keys`（打ち間違いの疑い）を持つ。
`message` と `details` にはファイル名・設定キー・環境変数名・設定ファイルのパスだけが入り、
ライブラリの置き場所は入らない。
```

`docs/workflow/library.md` の

```
前提: なし。解決順は `file_path` の明示 → `ion_mode` に対応する環境変数
（`MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`。研究室の参照ライブラリはリポジトリの外に
置いてここで指す）→ データディレクトリの `*_Loaded.msp2.dbs` → 設定済みの環境変数
（両方あれば `MSP_AMBIGUOUS`）→ データディレクトリの `*.msp`（複数あれば
```

を

```
前提: なし。解決順は `file_path` の明示 → `ion_mode` に対応する極性の設定
（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` → `lipidmix.local.toml` の
`[library] msp_positive` / `msp_negative`。`lipidmix/core/user_config.py` の
`get_setting` が引く。研究室の参照ライブラリはリポジトリの外に置いてここで指す）→
データディレクトリの `*_Loaded.msp2.dbs` → 設定済みの極性
（両方あれば `MSP_AMBIGUOUS`）→ データディレクトリの `*.msp`（複数あれば
```

にする。

`docs/cli.md` の

```
| `--ion-mode {positive,negative}` | 環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` のどちらを使うか |
```

を

```
| `--ion-mode {positive,negative}` | どちらの極性の設定（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`、無ければ `lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`）を使うか |
```

にする。

Run: `C:/Python314/python.exe -m pytest tests/test_workflow_docs.py -q`
Expected: PASS（`docs/workflow/library.md` が挙げた `lipidmix/core/user_config.py` と `get_setting` は実在する）

- [ ] **Step 7: CLAUDE.md を直す**

環境変数の段落のうち

```
クライアント＝Use-LLLM は `payload` を置く）/ `MSDIAL_MSP_POS` `MSDIAL_MSP_NEG`
（研究室の参照ライブラリ `.msp` の置き場所。`library_load(ion_mode=...)` が読む）/
`LIPIDMIX_LIBRARY_CACHE_DIR`（照合用 SQLite キャッシュの置き場所。既定 `data/.library-cache`）。
```

を

```
クライアント＝Use-LLLM は `payload` を置く）/
`LIPIDMIX_LIBRARY_CACHE_DIR`（照合用 SQLite キャッシュの置き場所。既定 `data/.library-cache`）/
`LIPIDMIX_CONFIG`（下記の設定ファイルの場所の上書き。worktree とテストが使う）。

**外部資産の場所は設定ファイルか環境変数で指す**（spec 2026-10-06）。リポジトリ直下の
`lipidmix.local.toml`（追跡外。雛形 `lipidmix.example.toml`）の `[msdial] exe` `lbm` と
`[library] msp_positive` `msp_negative`、または環境変数 `MSDIAL_EXE` `MSDIAL_LBM`
`MSDIAL_MSP_POS` `MSDIAL_MSP_NEG`（環境変数が優先）。読むのは `lipidmix/core/user_config.py`
の `get_setting` だけで、各リゾルバは `os.environ` を直に読まない。設定ファイルは呼ばれる
たびに読むので再起動は要らない。テストは conftest がこれらの環境変数を消し
`LIPIDMIX_CONFIG` を存在しないパスへ向けている（手元の設定ファイルを掴ませない）。
```

にする。

研究室ライブラリの段落の `本体はリポジトリの外に置いて環境変数で指し、` を `本体はリポジトリの外に置いて設定ファイルか環境変数で指し、` にする。

pipeline の決まりごとの

```
  - **`core/atomic_io.py` と `core/process_control.py` は stdlib だけの leaf**。
```

を

```
  - **`core/atomic_io.py` と `core/process_control.py` と `core/user_config.py` は stdlib だけの leaf**。
```

にする。

Run: `C:/Python314/python.exe -m pytest tests/test_readme_links.py -q`
Expected: PASS（CLAUDE.md に数量表現を足していないこと）

- [ ] **Step 8: 全テストを流してからコミット（バックグラウンド）**

```bash
git add lipidmix.example.toml tests/test_example_config.py scripts/new-worktree.sh README.md USAGE.md docs/output_format/library.md docs/workflow/library.md docs/cli.md CLAUDE.md
git commit -m "docs(config): 設定ファイルの雛形と導入手順を足し、.env.example を消す" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" > "$SCRATCH/commit-t5.log" 2>&1; echo exit=$?; grep -E "passed|failed|FAILED" "$SCRATCH/commit-t5.log" | tail -3
```

（`git rm .env.example` は Step 2 で index に載っている。）

- [ ] **Step 9: 追跡外の記録と流れ図を更新する**

1. main ツリーの `C:\Users\yuu18\Metabolomix_with_LLM\docs\task.md` の末尾に、日付見出し `## 2026-10-06 設定ファイル（user_config）の実装` を足し、コミット SHA・全テストの結果・残り（`.mcp.json` などの環境変数を設定ファイルへ移すかはユーザー判断、`--no-ff` マージはユーザー判断）を書く。既存の節は書き換えない。
2. 同じく `docs/HISTRY.md` の末尾に、実装で判明した事実（`tomllib` が BOM を拒否すること、`"C:\new\tool.exe"` が構文として通ること、`MsdialExeNotFoundError` に `CONFIG_INVALID` を載せた理由）を日付見出しで足す。
3. vault の流れ図 `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md` で、Console 経路（`console_plan` / `pipeline_run`）と `library_load` の前提条件として `MSDIAL_EXE` / `MSDIAL_MSP_*` を挙げている箇所を「`lipidmix.local.toml` の `[msdial] exe` / `[library] msp_*` か環境変数」に直し、`CONFIG_INVALID` で止まる分岐を足す。末尾の出典行の日付を 2026-10-06 にする。
```
