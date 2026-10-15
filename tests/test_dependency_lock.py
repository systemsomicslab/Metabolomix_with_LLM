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
