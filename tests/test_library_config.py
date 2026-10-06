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
