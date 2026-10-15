# Python の仮想環境化と 1 コマンドセットアップ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 配布先ユーザーが `setup.ps1` / `setup.sh` を 1 回実行するだけで、版を固定した `.venv`・MCP クライアント設定・`lipidmix.local.toml` 雛形が揃うようにし、開発側の道具も同じ `.venv` へ移す。

**Architecture:** `pyproject.toml`（`package = false`）＋`uv.lock` を版の正準とし、pip 経路用に `uv export` で書き出したロックを追跡する。入口スクリプトは uv があれば `uv sync --locked`、無ければ標準 venv+pip で `.venv` を作り、共通の後半（import 自己確認・設定生成）を `.venv` の Python で走る `scripts/setup_finish.py` に寄せる。`server_update` は同じロックを pip か `uv pip` で入れ直す。

**Tech Stack:** Python 3.14（`requires-python >=3.13`）、uv 0.11、PowerShell 7 / POSIX sh、pytest。

**Spec:** `docs/superpowers/specs/2026-10-15-python-venv-setup-design.md`

## Global Constraints

- 作業はすべて worktree `C:/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-python-venv`（ブランチ `feat/python-venv`）で行う。
- `[tool.uv] package = false` を必ず保つ（editable install が worktree に main のコードを掴ませるため）。
- `requires-python = ">=3.13"`、`.python-version` は `3.14`、`mcp<2` を保持。
- 追跡するロック: `uv.lock` `requirements.lock.txt` `requirements-dev.lock.txt`。`requirements.txt` / `requirements-dev.txt` はそれぞれ `-r requirements.lock.txt` / `-r requirements-dev.lock.txt` の 1 行。
- ルート直下の `.py` は `server.py` `check.py` だけ（`tests/test_package_layout.py`）。新しい Python は `scripts/` に置く。
- README.md / CLAUDE.md に数量表現を書かない（`tests/test_readme_links.py`）。
- 既存の設定ファイル（`.mcp.json` `lipidmix.local.toml`）を setup が**上書きしない**。
- Claude Desktop の設定ファイルは書き換えない（貼り付け断片の表示だけ）。
- R は 1 コマンドに含めない（`Rscript r/setup.R` のまま）。
- テストは `<worktree>/.venv/Scripts/python.exe -m pytest ...`（Task 1 で作る）。それ以前は `C:/Python314/python.exe -m pytest ...`。
- コミットは pre-commit が全テストを回すので**バックグラウンド実行**し、出力をファイルへ落として `FAILED` を引く。実行中はファイルを編集しない。
- コミットメッセージ末尾: `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>`

## Review Focus

1. **setup 後の作業ツリーが汚れる** — `setup` が生成した `.mcp.json` が追跡外として残ると、配布先の `server_update` が `dirty_worktree` で拒否し続ける。期待: setup 直後の `git status --porcelain` が空。→ Task 1 の gitignore テストに `.mcp.json` を足す。
2. **既存の `.mcp.json` を壊す** — 開発者や利用者が env（MSP パス等）を書き込んだ `.mcp.json` を setup が書き換えると復元できない。期待: バイト不変で、差だけを表示。→ Task 2 のテスト。
3. **ロックの書き出し忘れ** — `uv.lock` だけ更新して `requirements.lock.txt` を書き出し忘れると、uv 経路と pip 経路で版がずれる。期待: テストが落ちる。→ Task 1 のテスト。
4. **uv 製 venv に pip が無い** — `server_update` が `python -m pip` 固定のままだと uv 経路の利用者で必ず失敗する。期待: `uv pip` に切り替わり、どちらも無ければ `updated_with_warning` で手順を案内。→ Task 4 のテスト。
5. **Python 3.12 以下しか無い機械** — フォールバックが古い Python で venv を作ると import で謎の失敗をする。期待: 版を確かめて拒否し、uv か 3.13+ の導入を案内。→ Task 3 の手動検証（`-NoUv` で古い版を指せないため、版判定の式を Step で目視確認）。

---

### Task 1: 依存の宣言とロック

**Files:**
- Create: `pyproject.toml`
- Create: `.python-version`
- Create: `uv.lock`（`uv lock` が生成）
- Create: `requirements.lock.txt`, `requirements-dev.lock.txt`（`uv export` が生成）
- Modify: `requirements.txt`, `requirements-dev.txt`
- Modify: `.gitignore`
- Create: `tests/test_dependency_lock.py`
- Modify: `tests/test_gitignore_library.py`

**Interfaces:**
- Produces: `pyproject.toml` の `[project].dependencies`（Task 2 の import 自己確認が読む）、`requirements.lock.txt`（Task 3・4 が pip に渡す）、`<worktree>/.venv`（以降のテスト実行環境）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_dependency_lock.py`:

```python
"""依存の版の正準（uv.lock）と、pip 経路用の書き出しが食い違わないこと。

セットアップは uv があれば `uv sync --locked`、無ければ `pip install -r
requirements.lock.txt` で `.venv` を作る（spec 2026-10-15）。uv.lock だけ更新して
書き出しを忘れると、2 つの経路が別の版を入れる。uv が無くても走るよう、
uv.lock は TOML として直に読む。
"""
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_PIN = re.compile(r"^([A-Za-z0-9][A-Za-z0-9._-]*)==([^\s;\\]+)")


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def _locked_versions() -> dict[str, str]:
    data = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    return {_norm(p["name"]): p["version"] for p in data["package"] if "version" in p}


def _pins(filename: str) -> dict[str, str]:
    pins = {}
    for line in (ROOT / filename).read_text(encoding="utf-8").splitlines():
        match = _PIN.match(line.strip())
        if match:
            pins[_norm(match.group(1))] = match.group(2)
    return pins


def _direct_dependencies() -> list[str]:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return [_norm(re.split(r"[<>=!~;\[ ]", spec, maxsplit=1)[0])
            for spec in project["project"]["dependencies"]]


def test_runtime_lock_matches_uv_lock():
    locked = _locked_versions()
    pins = _pins("requirements.lock.txt")
    assert pins, "requirements.lock.txt に固定行が無い"
    mismatched = {n: (v, locked.get(n)) for n, v in pins.items() if locked.get(n) != v}
    assert mismatched == {}, "uv.lock と食い違う。`uv export` で書き出し直す"


def test_dev_lock_matches_uv_lock():
    locked = _locked_versions()
    pins = _pins("requirements-dev.lock.txt")
    mismatched = {n: (v, locked.get(n)) for n, v in pins.items() if locked.get(n) != v}
    assert mismatched == {}
    assert "pytest" in pins


def test_runtime_lock_covers_direct_dependencies_and_excludes_dev():
    pins = _pins("requirements.lock.txt")
    missing = [name for name in _direct_dependencies() if name not in pins]
    assert missing == []
    assert "pytest" not in pins


def test_requirements_files_are_pointers_to_the_locks():
    """旧版の server_update は requirements.txt を読む。中身は lock への委譲だけにする。"""
    assert (ROOT / "requirements.txt").read_text(encoding="utf-8").split() == [
        "-r", "requirements.lock.txt"]
    assert (ROOT / "requirements-dev.txt").read_text(encoding="utf-8").split() == [
        "-r", "requirements-dev.lock.txt"]


def test_project_itself_is_not_installed():
    """editable install は .pth で main ツリーを指し、worktree に main のコードを掴ませる。"""
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert project["tool"]["uv"]["package"] is False
```

`tests/test_gitignore_library.py` の末尾に追加:

```python
#: 追跡しなければならない（`*.txt` の無視に巻き込まれてはいけない）依存のロック。
_MUST_BE_TRACKABLE = ["requirements.txt", "requirements-dev.txt",
                      "requirements.lock.txt", "requirements-dev.lock.txt"]


@pytest.mark.skipif(not _git_available(), reason="git の作業ツリーでない（配備先など）")
@pytest.mark.parametrize("relpath", _MUST_BE_TRACKABLE)
def test_dependency_locks_are_not_ignored(relpath):
    result = subprocess.run(["git", "check-ignore", "--no-index", "-q", relpath], cwd=ROOT)
    assert result.returncode == 1, f"{relpath} が .gitignore で無視されている"


@pytest.mark.skipif(not _git_available(), reason="git の作業ツリーでない（配備先など）")
@pytest.mark.parametrize("relpath", [".mcp.json", ".venv/pyvenv.cfg"])
def test_setup_outputs_are_ignored(relpath):
    """setup が作るものが追跡外で残ると、server_update が dirty_worktree で拒否し続ける。"""
    result = subprocess.run(["git", "check-ignore", "--no-index", "-q", relpath], cwd=ROOT)
    assert result.returncode == 0, f"{relpath} が .gitignore で無視されていない"
```

- [ ] **Step 2: 失敗を確認**

Run: `C:/Python314/python.exe -m pytest tests/test_dependency_lock.py tests/test_gitignore_library.py -q`
Expected: `test_dependency_lock.py` は `FileNotFoundError`（uv.lock / pyproject.toml 不在）で FAIL、`requirements.lock.txt` の trackable 判定と `.mcp.json` の ignore 判定が FAIL。

- [ ] **Step 3: 宣言を書く**

`pyproject.toml`:

```toml
# 依存の宣言。版の正準は uv.lock（pip 経路用には requirements*.lock.txt へ書き出す）。
# 更新手順: 依存をここで直す → `uv lock` → 下の 2 行で書き出す → コミット。
#   uv export --frozen --no-hashes --no-emit-project --no-dev -o requirements.lock.txt
#   uv export --frozen --no-hashes --no-emit-project --group dev -o requirements-dev.lock.txt
[project]
name = "ms-data-parser"
version = "0.1.0"
description = "MCP server that parses and analyzes MS-DIAL outputs"
requires-python = ">=3.13"
dependencies = [
    "mcp<2",
    "msgpack",
    "lz4",
    "numpy",
    "matplotlib",
    "pandas",
    "scikit-learn",
    "scipy",
    "pygoslin",
    "httpx",
    "plotly",
    "rdkit",
    "PyMuPDF",
]

[dependency-groups]
dev = ["pytest"]

[tool.uv]
# プロジェクト自身は入れない。editable install の .pth が main ツリーを指すと、
# worktree から起動した server.py / pytest が main のコードを拾いうる
# （全 worktree が main ツリーの .venv を共有する前提。spec 2026-10-15）。
package = false
```

`.python-version`:

```
3.14
```

`requirements.txt`:

```
-r requirements.lock.txt
```

`requirements-dev.txt`:

```
-r requirements-dev.lock.txt
```

`.gitignore`: 先頭の `.venv/` の次の行に `.mcp.json` を、末尾の `*.txt` 行の直後のコメントを置き換えて例外を足す:

```gitignore
.venv/
# setup が生成する MCP クライアント設定（利用者ごとの絶対パスが入る）。
# 追跡外のまま残ると server_update が dirty_worktree で拒否する。
.mcp.json
```

```gitignore
*.txt
# 依存の宣言とロックは追跡する（*.txt の無視の例外）。
!requirements*.txt
```

- [ ] **Step 4: ロックを作り書き出す**

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-python-venv
uv lock
uv export --frozen --no-hashes --no-emit-project --no-dev -o requirements.lock.txt
uv export --frozen --no-hashes --no-emit-project --group dev -o requirements-dev.lock.txt
uv sync --locked --group dev
```

Expected: `.venv/Scripts/python.exe` ができる。`requirements.lock.txt` に `numpy==` 等の行があり `pytest` が無い。

- [ ] **Step 5: テストが通ることと、全スイートが .venv で緑であることを確認**

Run: `.venv/Scripts/python.exe -m pytest tests/test_dependency_lock.py tests/test_gitignore_library.py -q`
Expected: PASS

Run: `.venv/Scripts/python.exe -m pytest tests -q`
Expected: 全数 PASS（新しい版の依存でも既存テストが通る。落ちたら依存の版差が原因かを調べ、必要なら `pyproject.toml` に上限を足して `uv lock` からやり直す）。

- [ ] **Step 6: コミット（バックグラウンド）**

```bash
git add pyproject.toml .python-version uv.lock requirements.lock.txt requirements-dev.lock.txt requirements.txt requirements-dev.txt .gitignore tests/test_dependency_lock.py tests/test_gitignore_library.py
git commit -m "build: pyproject.toml と uv.lock で依存の版を固定し、pip 経路用のロックを書き出す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task1.log 2>&1
```

---

### Task 2: `scripts/setup_finish.py`（共通の後半）

**Files:**
- Create: `scripts/setup_finish.py`
- Create: `tests/test_setup_finish.py`

**Interfaces:**
- Consumes: `pyproject.toml` の `[project].dependencies`（Task 1）。
- Produces（Task 3 が呼ぶ）: `python scripts/setup_finish.py [--root <dir>]`。終了コード 0＝成功、1＝import 自己確認の失敗。モジュール関数:
  - `IMPORT_NAMES: dict[str, str]`（正規化した配布名 → import 名。載っていなければ配布名の `-` を `_` にしたもの）
  - `direct_dependencies(root: Path) -> list[str]`
  - `missing_imports(distributions: list[str], importer=importlib.import_module) -> list[str]`
  - `venv_python(root: Path) -> Path`
  - `server_entry(root: Path) -> dict`
  - `ensure_local_config(root: Path) -> str`（`"created"` / `"kept"` / `"no_example"`）
  - `ensure_mcp_json(root: Path) -> tuple[str, str]`（状態 `"created"` / `"kept"` / `"differs"`、利用者向けの文）
  - `desktop_snippet(root: Path) -> str`
  - `main(argv: list[str] | None = None) -> int`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_setup_finish.py`:

```python
"""setup.ps1 / setup.sh の共通の後半（scripts/setup_finish.py）を固定する。

縛るのは「利用者の手元を壊さない」こと: 既存の .mcp.json と lipidmix.local.toml は
バイト単位で触らず、違いがあれば言うだけにする（spec 2026-10-15 §3）。
"""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("setup_finish", ROOT / "scripts" / "setup_finish.py")
setup_finish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(setup_finish)


@pytest.fixture
def repo(tmp_path):
    (tmp_path / "pyproject.toml").write_text(
        '[project]\nname = "x"\ndependencies = ["mcp<2", "scikit-learn", "PyMuPDF"]\n',
        encoding="utf-8")
    (tmp_path / "lipidmix.example.toml").write_text("[msdial]\n", encoding="utf-8")
    (tmp_path / "server.py").write_text("", encoding="utf-8")
    return tmp_path


def test_direct_dependencies_strip_specifiers(repo):
    assert setup_finish.direct_dependencies(repo) == ["mcp", "scikit-learn", "pymupdf"]


def test_missing_imports_maps_distribution_to_module_names():
    seen = []

    def importer(name):
        seen.append(name)
        if name == "fitz":
            raise ImportError(name)

    missing = setup_finish.missing_imports(["scikit-learn", "pymupdf", "mcp"], importer)
    assert seen == ["sklearn", "fitz", "mcp"]
    assert missing == ["pymupdf"]


def test_venv_python_is_absolute_and_inside_dot_venv(repo):
    python = setup_finish.venv_python(repo)
    assert python.is_absolute()
    assert python.parent.parent == repo / ".venv"


def test_creates_local_config_from_example(repo):
    assert setup_finish.ensure_local_config(repo) == "created"
    assert (repo / "lipidmix.local.toml").read_text(encoding="utf-8") == "[msdial]\n"


def test_keeps_existing_local_config_byte_for_byte(repo):
    (repo / "lipidmix.local.toml").write_bytes(b"[library]\nmsp_positive = 'x'\n")
    assert setup_finish.ensure_local_config(repo) == "kept"
    assert (repo / "lipidmix.local.toml").read_bytes() == b"[library]\nmsp_positive = 'x'\n"


def test_creates_mcp_json_pointing_at_venv_python(repo):
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "created"
    entry = json.loads((repo / ".mcp.json").read_text(encoding="utf-8"))["mcpServers"]["ms-data-parser"]
    assert entry["command"] == setup_finish.venv_python(repo).as_posix()
    assert entry["args"] == ["-u", (repo / "server.py").as_posix()]
    assert entry["env"] == {"PYTHONIOENCODING": "utf-8"}


def test_keeps_matching_mcp_json(repo):
    setup_finish.ensure_mcp_json(repo)
    before = (repo / ".mcp.json").read_bytes()
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "kept"
    assert (repo / ".mcp.json").read_bytes() == before


def test_reports_but_never_rewrites_a_differing_mcp_json(repo):
    original = json.dumps({"mcpServers": {"ms-data-parser": {
        "command": "C:/Python314/python.exe", "args": ["server.py"],
        "env": {"MSDIAL_MSP_POS": "x"}}}}).encode()
    (repo / ".mcp.json").write_bytes(original)
    state, message = setup_finish.ensure_mcp_json(repo)
    assert state == "differs"
    assert "C:/Python314/python.exe" in message
    assert setup_finish.venv_python(repo).as_posix() in message
    assert (repo / ".mcp.json").read_bytes() == original


def test_mcp_json_without_our_server_is_reported_not_rewritten(repo):
    original = b'{"mcpServers": {"other": {"command": "node"}}}'
    (repo / ".mcp.json").write_bytes(original)
    state, _ = setup_finish.ensure_mcp_json(repo)
    assert state == "differs"
    assert (repo / ".mcp.json").read_bytes() == original


def test_desktop_snippet_is_valid_json_with_absolute_paths(repo):
    snippet = json.loads(setup_finish.desktop_snippet(repo))
    entry = snippet["mcpServers"]["ms-data-parser"]
    assert Path(entry["command"]).is_absolute()
    assert Path(entry["args"][-1]).is_absolute()


def test_main_fails_when_a_dependency_is_missing(repo, monkeypatch, capsys):
    monkeypatch.setattr(setup_finish, "missing_imports", lambda names, importer=None: ["rdkit"])
    assert setup_finish.main(["--root", str(repo)]) == 1
    assert "rdkit" in capsys.readouterr().out
    assert not (repo / ".mcp.json").exists()


def test_main_succeeds_and_writes_both_files(repo, monkeypatch):
    monkeypatch.setattr(setup_finish, "missing_imports", lambda names, importer=None: [])
    assert setup_finish.main(["--root", str(repo)]) == 0
    assert (repo / ".mcp.json").is_file()
    assert (repo / "lipidmix.local.toml").is_file()
```

- [ ] **Step 2: 失敗を確認**

Run: `.venv/Scripts/python.exe -m pytest tests/test_setup_finish.py -q`
Expected: FAIL（`scripts/setup_finish.py` が無い）

- [ ] **Step 3: 実装する**

`scripts/setup_finish.py`:

```python
"""setup.ps1 / setup.sh の共通の後半。`.venv` の Python で走る。

    .venv の python scripts/setup_finish.py [--root <リポジトリ>]

やること（spec 2026-10-15 §3）:
1. pyproject.toml の直接依存がすべて import できるかを確かめる（できなければ 1 で終わる）。
2. lipidmix.local.toml が無ければ雛形から複写する。
3. .mcp.json が無ければ .venv の Python を指して書く。あれば**触らず**、違いだけを言う。
4. Claude Desktop 用の貼り付け断片と、次の手順を表示する。

既存の設定ファイルを書き換えないのは、利用者が env（参照ライブラリのパス等）を
書き込んでいる可能性があり、上書きすると復元できないため。
"""
from __future__ import annotations

import argparse
import importlib
import json
import re
import shutil
import sys
import tomllib
from pathlib import Path

SERVER_NAME = "ms-data-parser"

#: 配布名（正規化済み）と import 名が違うものだけを載せる。
IMPORT_NAMES = {
    "scikit-learn": "sklearn",
    "pymupdf": "fitz",
}


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def direct_dependencies(root: Path) -> list[str]:
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    return [_norm(re.split(r"[<>=!~;\[ ]", spec, maxsplit=1)[0])
            for spec in project["project"]["dependencies"]]


def missing_imports(distributions: list[str], importer=importlib.import_module) -> list[str]:
    missing = []
    for dist in distributions:
        try:
            importer(IMPORT_NAMES.get(dist, dist.replace("-", "_")))
        except ImportError:
            missing.append(dist)
    return missing


def venv_python(root: Path) -> Path:
    venv = root.resolve() / ".venv"
    if sys.platform == "win32":
        return venv / "Scripts" / "python.exe"
    return venv / "bin" / "python"


def server_entry(root: Path) -> dict:
    return {
        "type": "stdio",
        "command": venv_python(root).as_posix(),
        "args": ["-u", (root.resolve() / "server.py").as_posix()],
        "env": {"PYTHONIOENCODING": "utf-8"},
    }


def ensure_local_config(root: Path) -> str:
    target = root / "lipidmix.local.toml"
    if target.exists():
        return "kept"
    example = root / "lipidmix.example.toml"
    if not example.is_file():
        return "no_example"
    shutil.copyfile(example, target)
    return "created"


def ensure_mcp_json(root: Path) -> tuple[str, str]:
    target = root / ".mcp.json"
    wanted = server_entry(root)
    if not target.exists():
        target.write_text(json.dumps({"mcpServers": {SERVER_NAME: wanted}}, indent=2) + "\n",
                          encoding="utf-8")
        return "created", f"{target} を作りました。"
    try:
        current = json.loads(target.read_text(encoding="utf-8"))["mcpServers"][SERVER_NAME]
    except (ValueError, KeyError, TypeError):
        return "differs", (f"{target} は既にありますが {SERVER_NAME} の項がありません（触っていません）。"
                           f"次の項を mcpServers に足してください:\n"
                           f"{json.dumps({SERVER_NAME: wanted}, indent=2)}")
    if current.get("command") == wanted["command"]:
        return "kept", f"{target} は既に .venv の Python を指しています。"
    return "differs", (f"{target} は既にあります（触っていません）。{SERVER_NAME} の command を\n"
                       f"  現在: {current.get('command')}\n"
                       f"  推奨: {wanted['command']}\n"
                       f"に手で直してください（env などほかの項はそのままで構いません）。")


def desktop_snippet(root: Path) -> str:
    return json.dumps({"mcpServers": {SERVER_NAME: server_entry(root)}}, indent=2)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    root = parser.parse_args(argv).root.resolve()

    missing = missing_imports(direct_dependencies(root))
    if missing:
        print("[setup] import できない依存があります: " + ", ".join(missing))
        print("[setup] setup をもう一度実行してください。直らなければ出力を添えて報告してください。")
        return 1
    print("[setup] 依存の import を確認しました。")

    config = ensure_local_config(root)
    if config == "created":
        print(f"[setup] {root / 'lipidmix.local.toml'} を作りました。MS-DIAL や参照ライブラリを"
              "使うときは、ここに場所を書いてください（任意）。")
    elif config == "kept":
        print("[setup] lipidmix.local.toml は既にあるので触っていません。")

    _, message = ensure_mcp_json(root)
    print(f"[setup] {message}")

    print("\n[setup] Claude Desktop で使うときは、claude_desktop_config.json の mcpServers に"
          "次を足してください\n（同じ名前のサーバを .mcp.json と二重に登録しないでください）:")
    print(desktop_snippet(root))
    print("\n[setup] 完了。MCP クライアントを再起動してください。")
    print("[setup] パスウェイ投影（pathway_*）を使うときだけ、R を入れて `Rscript r/setup.R` を 1 度実行します。")
    print("[setup] 更新は会話で server_update を呼ぶか、git pull の後に setup をもう一度実行します。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: テストが通ることを確認**

Run: `.venv/Scripts/python.exe -m pytest tests/test_setup_finish.py tests/test_package_layout.py -q`
Expected: PASS

Run（実物で自己確認が通るか）: `.venv/Scripts/python.exe scripts/setup_finish.py --root <scratchpad の空ディレクトリに pyproject.toml と server.py を複写したもの>`
Expected: 「依存の import を確認しました。」と終了コード 0。

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add scripts/setup_finish.py tests/test_setup_finish.py
git commit -m "feat(setup): import の自己確認と MCP 設定・雛形の生成を行う setup_finish.py

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task2.log 2>&1
```

---

### Task 3: 入口スクリプト `setup.ps1` / `setup.sh`

**Files:**
- Create: `setup.ps1`
- Create: `setup.sh`

**Interfaces:**
- Consumes: `uv.lock` `requirements*.lock.txt`（Task 1）、`scripts/setup_finish.py`（Task 2）。
- Produces（Task 6 の文書が案内する）: `./setup.ps1 [-Dev] [-NoUv]` / `sh setup.sh [--dev] [--no-uv]`。成功で 0、失敗で非 0 と 1 行の理由。

- [ ] **Step 1: `setup.ps1` を書く**

```powershell
<#
.SYNOPSIS
  ms-data-parser の Python 環境（.venv）と MCP 設定を 1 回で揃える。

.DESCRIPTION
  リポジトリ直下で実行する:  ./setup.ps1        （利用者）
                             ./setup.ps1 -Dev   （開発者。pytest も入る）
  uv があれば uv.lock どおりに入れる。無ければ Python 3.13 以上で標準の venv を作り、
  requirements.lock.txt（uv.lock の書き出し）を pip で入れる。最後に
  scripts/setup_finish.py が import の自己確認と .mcp.json / lipidmix.local.toml の生成を行う。
  R（パスウェイ投影）は含まない: 使うときだけ Rscript r/setup.R。

.PARAMETER Dev
  開発用の依存（pytest）も入れる。

.PARAMETER NoUv
  uv があっても使わず、標準の venv + pip で作る（検証用）。
#>
param([switch]$Dev, [switch]$NoUv)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

function Fail([string]$message) {
    Write-Host "[setup] 失敗: $message" -ForegroundColor Red
    exit 1
}

function Find-Python {
    # 3.13 以上の Python を探す。見つかれば起動に使う引数列を返す。
    $candidates = @(@('py', '-3.14'), @('py', '-3.13'), @('python'), @('python3'))
    foreach ($candidate in $candidates) {
        if (-not (Get-Command $candidate[0] -ErrorAction SilentlyContinue)) { continue }
        $exe = $candidate[0]
        $rest = @($candidate | Select-Object -Skip 1)
        & $exe @rest -c "import sys; sys.exit(0 if sys.version_info >= (3, 13) else 1)" 2>$null
        if ($LASTEXITCODE -eq 0) { return ,$candidate }
    }
    return $null
}

$uv = if ($NoUv) { $null } else { Get-Command uv -ErrorAction SilentlyContinue }
if ($uv) {
    Write-Host '[setup] uv で .venv を作ります（uv.lock どおり）。'
    $group = if ($Dev) { @('--group', 'dev') } else { @('--no-dev') }
    & uv sync --locked @group
    if ($LASTEXITCODE -ne 0) { Fail 'uv sync が失敗しました。上の出力を確認してください。' }
} else {
    $python = Find-Python
    if (-not $python) {
        Fail ('Python 3.13 以上も uv も見つかりません。どちらかを入れてから再実行してください。' +
              ' uv: https://docs.astral.sh/uv/getting-started/installation/ ' +
              ' Python: https://www.python.org/downloads/')
    }
    Write-Host "[setup] uv が無いので標準の venv で作ります（$($python -join ' ')）。"
    $exe = $python[0]
    $rest = @($python | Select-Object -Skip 1)
    & $exe @rest -m venv .venv
    if ($LASTEXITCODE -ne 0) { Fail 'python -m venv が失敗しました。' }
    $lock = if ($Dev) { 'requirements-dev.lock.txt' } else { 'requirements.lock.txt' }
    & .\.venv\Scripts\python.exe -m pip install --disable-pip-version-check -r $lock
    if ($LASTEXITCODE -ne 0) { Fail "pip install -r $lock が失敗しました。上の出力を確認してください。" }
}

& .\.venv\Scripts\python.exe scripts\setup_finish.py
if ($LASTEXITCODE -ne 0) { Fail 'セットアップの確認で問題が見つかりました（上の出力）。' }
```

- [ ] **Step 2: `setup.sh` を書く**

```sh
#!/bin/sh
# ms-data-parser の Python 環境（.venv）と MCP 設定を 1 回で揃える（mac / Linux / Git Bash）。
#
#   sh setup.sh          利用者
#   sh setup.sh --dev    開発者（pytest も入る）
#   sh setup.sh --no-uv  uv があっても標準の venv + pip で作る（検証用）
#
# uv があれば uv.lock どおりに入れる。無ければ Python 3.13 以上で venv を作り、
# requirements.lock.txt（uv.lock の書き出し）を pip で入れる。最後に
# scripts/setup_finish.py が import の自己確認と .mcp.json / lipidmix.local.toml の生成を行う。
# R（パスウェイ投影）は含まない: 使うときだけ Rscript r/setup.R。
set -e
cd "$(dirname "$0")"

dev=0
use_uv=1
for arg in "$@"; do
    case "$arg" in
        --dev) dev=1 ;;
        --no-uv) use_uv=0 ;;
        *) echo "[setup] 不明な引数: $arg" >&2; exit 2 ;;
    esac
done

fail() { echo "[setup] 失敗: $1" >&2; exit 1; }

venv_python() {
    if [ -x .venv/Scripts/python.exe ]; then echo .venv/Scripts/python.exe
    else echo .venv/bin/python; fi
}

if [ "$use_uv" = 1 ] && command -v uv >/dev/null 2>&1; then
    echo "[setup] uv で .venv を作ります（uv.lock どおり）。"
    if [ "$dev" = 1 ]; then group="--group dev"; else group="--no-dev"; fi
    # shellcheck disable=SC2086
    uv sync --locked $group || fail "uv sync が失敗しました。上の出力を確認してください。"
else
    python=""
    for candidate in python3.14 python3.13 python3 python; do
        if command -v "$candidate" >/dev/null 2>&1 &&
           "$candidate" -c "import sys; sys.exit(0 if sys.version_info >= (3, 13) else 1)" 2>/dev/null; then
            python="$candidate"; break
        fi
    done
    [ -n "$python" ] || fail "Python 3.13 以上も uv も見つかりません。どちらかを入れてから再実行してください。 uv: https://docs.astral.sh/uv/getting-started/installation/  Python: https://www.python.org/downloads/"
    echo "[setup] uv が無いので標準の venv で作ります（$python）。"
    "$python" -m venv .venv || fail "python -m venv が失敗しました。"
    if [ "$dev" = 1 ]; then lock=requirements-dev.lock.txt; else lock=requirements.lock.txt; fi
    "$(venv_python)" -m pip install --disable-pip-version-check -r "$lock" ||
        fail "pip install -r $lock が失敗しました。上の出力を確認してください。"
fi

"$(venv_python)" scripts/setup_finish.py || fail "セットアップの確認で問題が見つかりました（上の出力）。"
```

- [ ] **Step 3: uv 経路を新品の worktree で手動検証**

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM
git worktree add .worktrees/verify-uv feat/python-venv --detach
cd .worktrees/verify-uv
pwsh -File ./setup.ps1 -Dev
.venv/Scripts/python.exe -m pytest tests -q
git status --porcelain
```

Expected: setup が「依存の import を確認しました」「.mcp.json を作りました」で 0 終了。pytest 全数 PASS（pipeline / Console worker のテストが venv 下で通る＝`process_control` のリダイレクタ処理の裏取り）。`git status --porcelain` が空（Review Focus 1）。

- [ ] **Step 4: venv フォールバックを新品の worktree で手動検証**

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM
git worktree add .worktrees/verify-pip feat/python-venv --detach
cd .worktrees/verify-pip
pwsh -File ./setup.ps1 -NoUv
.venv/Scripts/python.exe -c "import sys; print(sys.version)"
git status --porcelain
```

Expected: 「uv が無いので標準の venv で作ります（py -3.14）」→ import 確認が通り 0 終了。`git status --porcelain` が空。
版判定の目視確認（Review Focus 5）: `py -3.14 -c "import sys; sys.exit(0 if sys.version_info >= (3, 13) else 1)"; echo $LASTEXITCODE` が 0、`(3, 99)` に書き換えた式で 1 になること。

- [ ] **Step 5: `setup.sh` を Git Bash で手動検証**

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM/.worktrees/verify-uv
rm -rf .venv .mcp.json lipidmix.local.toml
sh setup.sh
```

Expected: uv 経路で 0 終了、`.venv/Scripts/python.exe` を使って setup_finish が走る。mac / Linux 実機は無いので**未検証として HISTRY に記録**する（Task 6）。

- [ ] **Step 6: 検証用 worktree を片付け、コミット（バックグラウンド）**

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM
git worktree remove --force .worktrees/verify-uv
git worktree remove --force .worktrees/verify-pip
cd .worktrees/feat-python-venv
git add setup.ps1 setup.sh
git commit -m "feat(setup): uv を優先し無ければ標準 venv で .venv を作る setup.ps1 / setup.sh

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task3.log 2>&1
```

---

### Task 4: `server_update` を lock と uv 製 venv に対応させる

**Files:**
- Modify: `metabolomix/core/version.py`（定数 `REQUIREMENTS_FILE` 周辺、`_pip_install`、`apply_update` の依存判定と失敗文）
- Modify: `tests/test_server_update.py`

**Interfaces:**
- Consumes: `requirements.lock.txt`（Task 1）。
- Produces:
  - `DEPENDENCY_FILES: tuple[str, ...] = ("requirements.txt", "requirements.lock.txt", "uv.lock")`
  - `LOCK_FILE = "requirements.lock.txt"`
  - `_install_command(lock: Path) -> list[str] | None`
  - `_pip_install() -> tuple[bool, str]`（既存の名前と戻り値の形を保つ。既存テストの `no_pip` が差し替える）
  - `REQUIREMENTS_FILE` は削除（参照は version.py の中だけ）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_server_update.py` の「依存の更新」節の末尾に追加:

```python
@pytest.mark.parametrize("changed", ["uv.lock", "requirements.lock.txt"])
def test_installs_dependencies_when_lock_changed(monkeypatch, changed):
    """依存の正準は uv.lock。requirements.txt はロックへのポインタにすぎない。"""
    monkeypatch.setattr(version, "_git", FakeGit(_with(
        **{"diff__--name-only": f"{changed}\nserver.py"})))
    calls = []
    monkeypatch.setattr(version, "_pip_install",
                        lambda: (calls.append(True), (True, ""))[1])

    result = version.apply_update()

    assert result["dependencies"] == "installed"
    assert len(calls) == 1


def test_install_command_uses_pip_when_available(monkeypatch, tmp_path):
    monkeypatch.setattr(version.importlib.util, "find_spec",
                        lambda name: object() if name == "pip" else None)
    lock = tmp_path / "requirements.lock.txt"
    assert version._install_command(lock) == [
        version.sys.executable, "-m", "pip", "install", "-r", str(lock)]


def test_install_command_falls_back_to_uv_pip_in_a_uv_venv(monkeypatch, tmp_path):
    """uv が作った venv は pip を持たない。"""
    monkeypatch.setattr(version.importlib.util, "find_spec", lambda name: None)
    monkeypatch.setattr(version.shutil, "which", lambda name: "C:/bin/uv.exe" if name == "uv" else None)
    lock = tmp_path / "requirements.lock.txt"
    assert version._install_command(lock) == [
        "C:/bin/uv.exe", "pip", "install", "--python", version.sys.executable, "-r", str(lock)]


def test_pip_install_reports_failure_when_no_installer(monkeypatch, tmp_path):
    (tmp_path / "requirements.lock.txt").write_text("numpy==2.0.0\n", encoding="utf-8")
    monkeypatch.setattr(version, "_repo_root", lambda: tmp_path)
    monkeypatch.setattr(version, "_install_command", lambda lock: None)

    ok, log = version._pip_install()

    assert ok is False
    assert "uv" in log and "pip" in log


def test_dependency_failure_points_to_setup_script(monkeypatch):
    monkeypatch.setattr(version, "_git", FakeGit(_with(
        **{"diff__--name-only": "uv.lock"})))
    monkeypatch.setattr(version, "_pip_install", lambda: (False, "boom"))

    result = version.apply_update()

    assert result["status"] == "updated_with_warning"
    assert "setup.ps1" in result["message"]
```

- [ ] **Step 2: 失敗を確認**

Run: `.venv/Scripts/python.exe -m pytest tests/test_server_update.py -q`
Expected: 新しい 6 件が FAIL（`uv.lock` で入れ直さない／`_install_command` 不在／`importlib` `shutil` 未 import／文言に `setup.ps1` が無い）。既存は PASS。

- [ ] **Step 3: 実装する**

`metabolomix/core/version.py`:

import 節に `import importlib.util` と `import shutil` を足す。

定数の置き換え（`REQUIREMENTS_FILE = "requirements.txt"` とその上のコメント行を削除して）:

```python
#: 依存の宣言とロック。pull でどれかが変わったときだけ入れ直す。正準は uv.lock で、
#: requirements.txt はロックへのポインタ（旧版の server_update が読むので残っている）。
DEPENDENCY_FILES = ("requirements.txt", "requirements.lock.txt", "uv.lock")
#: 入れ直しに使うロック（uv.lock の pip 向け書き出し。spec 2026-10-15）。
LOCK_FILE = "requirements.lock.txt"
#: 依存の入れ直しに失敗したときに案内する手順。
SETUP_COMMAND = "./setup.ps1（mac / Linux は sh setup.sh）"
```

`_pip_install` を置き換え、`_install_command` を足す:

```python
def _install_command(lock: Path) -> list[str] | None:
    """`lock` を走行中の Python へ入れるコマンド。入れる手段が無ければ None。

    `sys.executable` を使うのは、PATH 上の `pip` が別の Python を指しうるため。
    uv が作った `.venv` は pip を持たないので、そのときは `uv pip --python` で
    同じ解釈系へ入れる。`uv sync` は使わない（dev group の pytest など、
    ロックに無いものを消してしまうため）。
    """
    if importlib.util.find_spec("pip") is not None:
        return [sys.executable, "-m", "pip", "install", "-r", str(lock)]
    uv = shutil.which("uv")
    if uv:
        return [uv, "pip", "install", "--python", sys.executable, "-r", str(lock)]
    return None


def _pip_install() -> tuple[bool, str]:
    """ロックを現在の Python へ入れ直す。`(成功したか, 出力)`。"""
    lock = _repo_root() / LOCK_FILE
    if not lock.is_file():
        return True, ""
    command = _install_command(lock)
    if command is None:
        return False, "pip も uv も見つからないため依存を入れられませんでした。"
    try:
        out = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL, timeout=PIP_TIMEOUT_SEC,
            text=True, errors="replace")
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"{type(exc).__name__}: {exc}"
    return out.returncode == 0, (out.stdout or "").strip()
```

`apply_update` の判定:

```python
    requirements_changed = any(
        line.strip() in DEPENDENCY_FILES for line in changed.splitlines())
```

失敗文の `手元で \`pip install -r {REQUIREMENTS_FILE}\` を実行してください。` を
`手元のリポジトリ直下で \`{SETUP_COMMAND}\` を実行してください。` に置き換える。

- [ ] **Step 4: テストが通ることを確認**

Run: `.venv/Scripts/python.exe -m pytest tests/test_server_update.py tests/test_version_update.py -q`
Expected: PASS（既存の `requirements.txt` 変化のテストも `DEPENDENCY_FILES` に含まれるので通る）

Run: `grep -rn REQUIREMENTS_FILE metabolomix tests`
Expected: 該当なし

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add metabolomix/core/version.py tests/test_server_update.py
git commit -m "fix(update): server_update が uv.lock の変化で入れ直し、pip の無い uv 製 venv では uv pip を使う

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task4.log 2>&1
```

---

### Task 5: 開発側の解釈系を `.venv` へ切り替える

**Files:**
- Modify: `.githooks/pre-commit`
- Modify: `scripts/new-worktree.sh`
- Modify: `.vscode/mcp.json`（追跡済み。`.gitignore` の `.vscode/` は追跡済みファイルに効かない）

**Interfaces:**
- Consumes: `.venv` の配置（Task 1・3）。
- Produces: 開発者は main ツリーで 1 度 `./setup.ps1 -Dev` を実行すれば、pre-commit と新しい worktree の `.mcp.json` がその `.venv` を使う。

- [ ] **Step 1: pre-commit の解釈系選択を書き換える**

`.githooks/pre-commit` の

```sh
PY=C:/Python314/python.exe
[ -f "$PY" ] || PY=python
```

を次に置き換える:

```sh
# 解釈系: 自ツリーの .venv → main ツリーの .venv（worktree は自前の .venv を持たない）
# → 従来の C:/Python314 → PATH の python（spec 2026-10-15）。
main=$(cd "$(dirname "$(git rev-parse --git-common-dir)")" && pwd)
PY=""
for candidate in "$root/.venv/Scripts/python.exe" "$root/.venv/bin/python" \
                 "$main/.venv/Scripts/python.exe" "$main/.venv/bin/python" \
                 C:/Python314/python.exe; do
    if [ -f "$candidate" ]; then PY=$candidate; break; fi
done
[ -n "$PY" ] || PY=python
echo "pre-commit: interpreter $PY"
```

- [ ] **Step 2: new-worktree.sh の command を main の .venv にする**

`scripts/new-worktree.sh` の `config_env=""` の直前に追加:

```sh
# MCP サーバは main ツリーの .venv で動かす（pyproject の package = false により
# import は常にこの worktree の server.py の隣から解決されるので共有して安全）。
python_w="C:/Python314/python.exe"
for candidate in "$main/.venv/Scripts/python.exe" "$main/.venv/bin/python"; do
    if [ -f "$candidate" ]; then python_w=$(win "$candidate"); break; fi
done
```

ヒアドキュメント内の `"command": "C:/Python314/python.exe",` を `"command": "$python_w",` に置き換える。

- [ ] **Step 3: `.vscode/mcp.json` を .venv にする**

```json
{
	"servers": {
		"ms-data-parser": {
			"type": "stdio",
			"command": "${workspaceFolder}\\.venv\\Scripts\\python.exe",
			"args": [
				"-u",
				"${workspaceFolder}\\server.py"
			]
		}
	},
	"inputs": []
}
```

- [ ] **Step 4: 動作確認**

フックを丸ごと走らせると全テストが回るので、選択ロジックだけを同じ式で確かめる（実フックでの確認は Step 5 のコミットログ先頭行）:

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM/.worktrees/feat-python-venv
sh -c 'root=$(git rev-parse --show-toplevel); main=$(cd "$(dirname "$(git rev-parse --git-common-dir)")" && pwd); for c in "$root/.venv/Scripts/python.exe" "$main/.venv/Scripts/python.exe" C:/Python314/python.exe; do [ -f "$c" ] && { echo "$c"; break; }; done'
```

Expected: 自ツリーの `.venv/Scripts/python.exe`。

`new-worktree.sh` は main ツリーにまだ `.venv` が無いので、Task 7 で main に setup を走らせた後に確認する。

- [ ] **Step 5: コミット（バックグラウンド）**

```bash
git add .githooks/pre-commit scripts/new-worktree.sh .vscode/mcp.json
git commit -m "build: pre-commit・worktree・VS Code の MCP 設定を .venv の Python へ切り替える

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task5.log 2>&1
```

Expected: ログ先頭に `pre-commit: interpreter .../feat-python-venv/.venv/Scripts/python.exe`、全テスト PASS。

---

### Task 6: 文書の更新

**Files:**
- Modify: `README.md`（Quick start と開発者向けテスト節）
- Modify: `DEPLOY.md`（Member Setup / MCP Client Example / Receiving Updates）
- Modify: `CLAUDE.md`（コマンド節、Python の段落、並列節の `.mcp.json` の説明、踏みやすい罠の `*.txt`）
- Modify: `USAGE.md` §16（Python 側の手順に触れている箇所があれば `setup` へ寄せる。R の手順は変えない）

**Interfaces:**
- Consumes: Task 3 のコマンド名とフラグ、Task 4 の `SETUP_COMMAND` の文言。

- [ ] **Step 1: README の Quick start を書き換える**

`## Quick start` から `Point it at your data` の直前までを次に置き換える:

````markdown
## Quick start

From the repository root, run the setup script once. It uses [uv](https://docs.astral.sh/uv/)
when available (uv also fetches Python if you don't have it) and otherwise falls back to a
standard `venv` with Python 3.13+. Dependency versions are pinned by `uv.lock`.

```powershell
./setup.ps1          # Windows (PowerShell)
```

```bash
sh setup.sh          # macOS / Linux
```

It creates `.venv/`, checks that every dependency imports, copies `lipidmix.local.toml` from
the example if missing, writes a `.mcp.json` for Claude Code, and prints the snippet to paste
into Claude Desktop's `claude_desktop_config.json`. Existing `.mcp.json` and
`lipidmix.local.toml` are never overwritten — differences are reported instead.
````

開発者向けの節（`Tests run from the repository root.` 以降のコードブロック）を:

````markdown
Tests run from the repository root. They need `pytest`, which is not part of the
runtime dependencies; set up with the development group first:

```powershell
./setup.ps1 -Dev
.venv/Scripts/python.exe -m pytest tests -q
```
````

- [ ] **Step 2: DEPLOY.md を書き換える**

`## Member Setup` の最初のコードブロックまでを:

````markdown
## Member Setup

Run the setup script once from the clone (see [README.md](README.md) "Quick start"):

```powershell
./setup.ps1
```

It builds `.venv/` (uv if installed, otherwise Python 3.13+ `venv` + pip with
`requirements.lock.txt`), writes `.mcp.json`, and prints the Claude Desktop snippet.
````

その後の手動起動例の `python server.py` を `.venv\Scripts\python.exe server.py` に、`## MCP Client Example` の `"command": "python",` を `"command": "C:\\Users\\<name>\\Lipidmix_with_LLM\\.venv\\Scripts\\python.exe",` に置き換える。

`## Receiving Updates` の該当文を置き換える:
- `only when \`requirements.txt\` changed in that range, reinstalls dependencies with the same Python that runs the server.` → `only when \`uv.lock\` or \`requirements.lock.txt\` changed in that range, reinstalls the locked dependencies into the same Python that runs the server (with pip, or with \`uv pip\` in a uv-made \`.venv\`, which has no pip).`
- `run \`pip install -r requirements.txt\` by hand **before** restarting` → `re-run \`./setup.ps1\` (or \`sh setup.sh\`) **before** restarting`
- 手作業の等価コードブロックの 2 行目を `./setup.ps1` に。

`## Receiving Updates` の末尾（`The check stays silent ...` 段落の後）に追加:

```markdown
### Moving an existing clone to `.venv`

Clones set up before `.venv` existed run the server with a global Python. To move:
run `./setup.ps1` once, then change `command` in your MCP client configuration
(Claude Desktop's `claude_desktop_config.json`, or `.mcp.json`) to the `.venv` Python
the script prints. Keep any `env` entries you added. Until you switch, `server_update`
keeps working against the global Python.
```

- [ ] **Step 3: CLAUDE.md を書き換える**

コマンドのコードブロック:

```bash
./setup.ps1 -Dev                                   # 初回のみ（pwsh）。.venv を作り pytest も入れる
.venv/Scripts/python.exe -m pytest tests -q
```

`MCP サーバ起動: \`C:/Python314/python.exe server.py\`` → `MCP サーバ起動: \`.venv/Scripts/python.exe server.py\``。

`**Python は \`C:/Python314/python.exe\` を使う**` で始まる段落を丸ごと次に置き換える:

```markdown
**Python は `.venv` の Python を使う**（`./setup.ps1 -Dev` が作る。uv があれば `uv.lock` どおり、
無ければ標準 venv + `requirements-dev.lock.txt`。spec 2026-10-15）。依存を変えるときは
`pyproject.toml` を直し `uv lock` → `pyproject.toml` 冒頭の 2 行の `uv export` → コミット
（`tests/test_dependency_lock.py` が書き出し忘れを捕まえる）。`[tool.uv] package = false` は外さない
（editable install が worktree に main のコードを掴ませる）。worktree は自前の `.venv` を持たず
main ツリーの `.venv` を共有する（pre-commit と `new-worktree.sh` がそれを選ぶ）。
Claude Code 用の `.mcp.json` は追跡外（`.gitignore`）で、setup が無ければ作る。
追跡されている `.vscode/mcp.json` も `.venv` を指す。
```

並列節の「スクリプトが worktree 専用の `.mcp.json` を生成する。」の文の後に「`command` は main ツリーの `.venv` の Python。」を足す。

踏みやすい罠の `*.txt` 行を:

```markdown
- `.gitignore` に `*.txt` があるため、**新規の .txt は無言で追跡漏れする**（`git add -f` が要る）。
  例外は `requirements*.txt`（依存のロック）だけ。
```

- [ ] **Step 4: USAGE.md §16 を確認する**

Run: `grep -n "pip\|python\|requirements" USAGE.md`
§16 に Python の導入手順があれば「`./setup.ps1` を済ませた上で」に寄せる。R の手順（`Rscript r/setup.R`）は変えない。該当が無ければ変更しない。

- [ ] **Step 5: 腐敗防止テストを通す**

Run: `.venv/Scripts/python.exe -m pytest tests/test_readme_links.py tests/test_workflow_docs.py -q`
Expected: PASS（README / CLAUDE.md に数量表現を書いていないこと、リンクが実在すること）

- [ ] **Step 6: コミット（バックグラウンド）**

```bash
git add README.md DEPLOY.md CLAUDE.md USAGE.md
git commit -m "docs: セットアップを setup.ps1 / setup.sh と .venv に寄せ、既存クローンの移行手順を足す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task6.log 2>&1
```

---

### Task 7: 統合・実機確認・記録

**Files:**
- 追跡外: main ツリーの `.venv/`（新規）、`.venv-1/`（削除）、`.mcp.json`（`command` の手直し）
- 追跡外: main ツリーの `docs/HISTRY.md`（追記）、`docs/task.md`（項目の削除）
- リポジトリ外: `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`

- [ ] **Step 1: 全体レビューとマージ**

ブランチ全体のレビュー後、main ツリーで:

```bash
cd C:/Users/yuu18/Metabolomix_with_LLM
git checkout main
git merge --no-ff feat/python-venv -m "Merge feat/python-venv: Python の仮想環境化と 1 コマンドセットアップ

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/merge.log 2>&1
```

（main ツリーが別ブランチ `feat/mslm-remote-api` 上にあるので、checkout の前に作業ツリーがクリーンであることと、並列エージェントがいないことを確認する。いればユーザーに確認する。）

- [ ] **Step 2: main ツリーに .venv を作り、自分の .mcp.json を切り替える**

```powershell
cd C:/Users/yuu18/Metabolomix_with_LLM
./setup.ps1 -Dev
```

Expected: `.mcp.json` は既にあるので `differs` で「現在: C:/Python314/python.exe / 推奨: .../.venv/Scripts/python.exe」と表示される。`.mcp.json` の `command` だけを推奨値に手で直す（env は残す）。

- [ ] **Step 3: MCP 起動を確認**

Claude Code の `/mcp` で ms-data-parser を再接続し、`list_data_files` が応答すること。

- [ ] **Step 4: new-worktree.sh の確認**

```bash
sh scripts/new-worktree.sh tmp/verify-wt
grep command .worktrees/tmp-verify-wt/.mcp.json
git worktree remove .worktrees/tmp-verify-wt
git branch -D tmp/verify-wt
```

Expected: `command` が main ツリーの `.venv/Scripts/python.exe`。

- [ ] **Step 5: `.venv-1/` を削除**

```bash
rm -rf C:/Users/yuu18/Metabolomix_with_LLM/.venv-1
```

- [ ] **Step 6: 記録**

- `docs/HISTRY.md` 末尾に `## 2026-10-15 Python の仮想環境化（feat/python-venv）` を追記: 合意事項（spec の §1）、実施内容、検証結果（uv 経路・`-NoUv` 経路・Git Bash `setup.sh`・MCP 接続）、**mac / Linux 実機の `setup.sh` は未検証**、`.mcp.json` を `.gitignore` に入れた理由（dirty_worktree）。
- `docs/task.md` から「TODO（セットアップ・別 spec）: Python の仮想環境化」の行を削除。「起票者の PC では MCP サーバが `.venv\Scripts\python.exe` で起動している。CLAUDE.md の「`.venv/` は存在しない」…」の行も、記述を削除したので消す。（並列エージェントがいるなら末尾追記だけにして整理は後で。）
- vault の流れ図: セットアップの入口（利用者が最初に何をするか）を `setup.ps1` / `setup.sh` → `.venv` → MCP 設定 → （任意）`Rscript r/setup.R` に書き直し、末尾の出典行の日付を 2026-10-15 にする。

- [ ] **Step 7: ブランチと worktree の片付け**

```bash
git worktree remove .worktrees/feat-python-venv
git branch -d feat/python-venv
```

（`origin` 側は残す。push はユーザーに確認してから。）
