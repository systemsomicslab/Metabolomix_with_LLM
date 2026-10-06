"""外部資産の場所の設定（`metabolomix.core.user_config`）。spec 2026-10-06。

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

from metabolomix.core import user_config
from metabolomix.core.user_config import ConfigInvalidError

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
    from metabolomix.core import mcp_core
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
