"""研究室の参照ライブラリの指定（`library_load` の `file_path`）。

B-2b でライブラリの探索（データフォルダの `.msp` / `.dbs`・設定の極性の MSP）は消えた。
指定先は `file_path` の明示か、データセットが覚えたもの（`tests/test_pathless_tools.py`）。
置き場所（ディレクトリ）は戻り値に出さない規則を、明示のパスでも守る。
"""
from __future__ import annotations

import json

import pytest

from metabolomix.core.path_resolvers import LibraryPathError, resolve_library_path


@pytest.fixture
def lab(tmp_path):
    return tmp_path / "secret-lab-share"


def _msp(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("NAME: x\nPRECURSORMZ: 100\nNum Peaks: 0\n\n", encoding="utf-8")
    return path


def test_resolve_library_path_does_not_search(tmp_path, monkeypatch):
    """パスを省いても、データフォルダ・設定の極性の MSP を探さない（None を返すだけ）。"""
    from tests.config_fixtures import use_config

    pos = _msp(tmp_path / "lab" / "pos.msp")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(pos)})
    monkeypatch.chdir(_msp(tmp_path / "cwd" / "x.msp").parent)
    assert resolve_library_path() is None
    assert resolve_library_path(str(pos)) == str(pos)


def test_missing_explicit_file_names_the_file_not_the_directory(lab):
    with pytest.raises(LibraryPathError) as exc:
        resolve_library_path(str(lab / "neg.msp"))
    assert exc.value.code == "LIBRARY_NOT_FOUND"
    assert "neg.msp" in exc.value.message
    dumped = exc.value.message + json.dumps(exc.value.details, ensure_ascii=False)
    assert "secret-lab-share" not in dumped


def test_library_load_missing_file_does_not_leak_the_directory(lab):
    from metabolomix.library.tools import library_load

    raw = library_load(file_path=str(lab / "pos.msp"))
    payload = json.loads(raw)
    assert payload["code"] == "LIBRARY_NOT_FOUND"
    assert "pos.msp" in payload["message"]
    assert "secret-lab-share" not in raw


def test_library_load_open_failure_does_not_leak_a_forward_slash_directory(
        tmp_path, lab, monkeypatch):
    """`/` 区切りで書かれた置き場所は、OSError の文言にも `/` のまま出うる。伏せ漏れしない。"""
    from metabolomix.library import tools as library_tools
    from tests.config_fixtures import use_config

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    pos = _msp(lab / "pos.msp")

    def boom(path, **_kwargs):
        raise OSError(f"cannot open {lab.as_posix()}/pos.msp.sqlite")

    monkeypatch.setattr(library_tools, "open_store", boom)
    raw = library_tools.library_load(file_path=f"{lab.as_posix()}/pos.msp")
    assert json.loads(raw)["status"] == "error"
    assert pos.exists()
    assert lab.as_posix() not in raw
    assert "secret-lab-share" not in raw
