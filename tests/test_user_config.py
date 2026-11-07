"""外部資産の場所の設定（`metabolomix.core.user_config`）。spec 2026-10-06・2026-11-06 §7.4。

設定ファイル（`LIPIDMIX_CONFIG` → `<repo>/lipidmix.local.toml`）だけを読む。資産の場所を
上書きする環境変数は読まない。conftest が `LIPIDMIX_CONFIG` を存在しないパスに向けて
いるので、各テストは必要なものだけを tmp の設定ファイルに書く。
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


# ---------- 設定のキーの表 ----------

def test_settings_table_has_a_kind_per_key():
    assert user_config.SETTINGS == {
        "msdial.exe": "path",
        "library.msp_positive": "path",
        "library.msp_negative": "path",
        "library.lbm": "path",
        "mslipidmapper.rscript": "path",
        "mslipidmapper.browser": "path",
        "mslipidmapper.mslipidmapper_api_url": "text",
        "mslipidmapper.api_token": "text",
    }
    assert user_config.MOVED_KEYS == {"msdial.lbm": "library.lbm"}
    assert user_config.NON_PATH_SETTINGS == frozenset(
        {"mslipidmapper.mslipidmapper_api_url", "mslipidmapper.api_token"})


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


# ---------- 値の形 ----------

def test_reads_a_value_from_the_config_file(cfg, tmp_path):
    exe = tmp_path / "MSDIALCUI.exe"
    _write(cfg, f"[msdial]\nexe = '{exe}'\n")
    setting = user_config.get_setting("msdial.exe")
    assert setting == user_config.Setting(key="msdial.exe", value=str(exe), config_file=str(cfg))


@pytest.mark.parametrize("env_var", [
    "MSDIAL_EXE", "MSDIAL_LBM", "MSDIAL_MSP_POS", "MSDIAL_MSP_NEG", "MSDIAL_LIBRARY_LBM",
    "MSLIPIDMAPPER_RSCRIPT", "MSLIPIDMAPPER_BROWSER", "MSLIPIDMAPPER_API_URL",
    "MSLIPIDMAPPER_API_TOKEN"])
def test_env_var_is_ignored(cfg, tmp_path, monkeypatch, env_var):
    """spec §7.4: 資産の場所の環境変数は読まない。設定ファイルの値か None。"""
    monkeypatch.setenv(env_var, str(tmp_path / "from_env.exe"))
    for key in user_config.SETTINGS:
        assert user_config.get_setting(key) is None
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'from_file.exe'}'\n"
                f"[library]\nlbm = '{tmp_path / 'from_file.lbm2'}'\n"
                f"msp_positive = '{tmp_path / 'from_file.msp'}'\n"
                "[mslipidmapper]\nmslipidmapper_api_url = 'http://file'\n")
    assert user_config.get_setting("msdial.exe").value == str(tmp_path / "from_file.exe")
    assert user_config.get_setting("library.lbm").value == str(tmp_path / "from_file.lbm2")
    assert user_config.get_setting("library.msp_positive").value == str(tmp_path / "from_file.msp")
    assert user_config.get_setting("mslipidmapper.mslipidmapper_api_url").value == "http://file"
    assert user_config.get_setting("library.msp_negative") is None


def test_relative_path_resolves_against_the_config_folder(cfg):
    _write(cfg, "[library]\nmsp_positive = 'libs/pos.msp'\n")
    setting = user_config.get_setting("library.msp_positive")
    assert setting.value == str(cfg.parent / "libs" / "pos.msp")


def test_text_setting_is_stripped_and_not_resolved(cfg):
    _write(cfg, "[mslipidmapper]\nmslipidmapper_api_url = '  http://10.0.0.1:8000  '\n")
    assert user_config.get_setting("mslipidmapper.mslipidmapper_api_url").value == \
        "http://10.0.0.1:8000"


def test_old_msdial_lbm_key_is_rejected(cfg, tmp_path):
    """決めたこと 4: `[msdial] lbm` は `[library] lbm` へ移した。黙って無視すると Console の
    LBM の解決が変わる（設定の段を飛ばして実行体のフォルダへ落ちる）ので CONFIG_INVALID。"""
    _write(cfg, f"[msdial]\nlbm = '{tmp_path / 'auto.lbm2'}'\n")
    for key in ("library.lbm", "msdial.exe"):
        with pytest.raises(ConfigInvalidError) as exc:
            user_config.get_setting(key)
        assert exc.value.code == "CONFIG_INVALID"
        assert "[msdial] lbm は [library] lbm へ移しました" in exc.value.message
        assert exc.value.details()["config_file"] == str(cfg)
    assert "msdial.lbm" not in user_config.SETTINGS


def test_library_lbm_is_the_single_lbm_key(cfg, tmp_path):
    _write(cfg, f"[library]\nlbm = '{tmp_path / 'lab_lipids.lbm2'}'\n")
    setting = user_config.get_setting("library.lbm")
    assert setting.value == str(tmp_path / "lab_lipids.lbm2")
    assert setting.config_file == str(cfg)


def test_int_setting_rejects_bool_and_string(cfg, monkeypatch):
    """種類 int は TOML の整数だけ（bool は int の部分型だが拒む）。Task 9 の server.* が使う。"""
    monkeypatch.setitem(user_config.SETTINGS, "server.max_dataset_handles", "int")
    _write(cfg, "[server]\nmax_dataset_handles = 12\n")
    setting = user_config.get_setting("server.max_dataset_handles")
    assert setting.value == 12 and type(setting.value) is int
    for bad in ("true", "'12'", "1.5"):
        user_config._cache.clear()
        _write(cfg, f"[server]\nmax_dataset_handles = {bad}\n")
        with pytest.raises(ConfigInvalidError) as exc:
            user_config.get_setting("server.max_dataset_handles")
        assert "整数" in exc.value.message


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
    cfg.write_text(f"﻿[msdial]\nexe = '{tmp_path / 'a.exe'}'\n", encoding="utf-8")
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


def test_broken_config_file_is_config_invalid_even_with_an_env_var(cfg, monkeypatch):
    """環境変数で設定ファイルを迂回する道は無くなった。"""
    _write(cfg, "[msdial\n")
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    with pytest.raises(ConfigInvalidError):
        user_config.get_setting("msdial.exe")


# ---------- エラー文の材料 ----------

def test_describe_missing_without_a_config_file(cfg):
    assert user_config.describe_missing("msdial.exe") == {
        "setting": "msdial.exe", "config_key": "[msdial] exe", "source": None,
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


def test_describe_missing_reports_the_moved_key_without_raising(cfg, tmp_path):
    _write(cfg, f"[msdial]\nlbm = '{tmp_path / 'auto.lbm2'}'\n")
    assert "[library] lbm" in user_config.describe_missing("library.lbm")["config_invalid"]


def test_missing_hint_when_the_config_file_is_absent(cfg):
    hint = user_config.missing_hint("msdial.exe", "MS-DIAL Console の実行体")
    assert "lipidmix.example.toml" in hint
    assert "[msdial] exe" in hint
    assert "MSDIAL_EXE" not in hint
    assert "環境変数" not in hint


def test_missing_hint_when_the_config_file_exists(cfg):
    _write(cfg, "[library]\n")
    hint = user_config.missing_hint("msdial.exe", "MS-DIAL Console の実行体")
    assert hint.startswith("lipidmix.local.toml の [msdial] exe に")
    assert "環境変数" not in hint


def test_setting_label_names_the_config_file_and_key(cfg, tmp_path):
    _write(cfg, f"[msdial]\nexe = '{tmp_path / 'a.exe'}'\n")
    assert user_config.setting_label(user_config.get_setting("msdial.exe")) == \
        "lipidmix.local.toml の [msdial] exe"


# ---------- テスト用の書き手 ----------

def test_use_config_writes_merges_and_points_at_the_file(tmp_path, monkeypatch):
    from tests.config_fixtures import use_config
    exe = tmp_path / "app" / "MSDIALCUI.exe"
    path = use_config(monkeypatch, tmp_path, {"msdial.exe": exe})
    assert os.environ["LIPIDMIX_CONFIG"] == str(path)
    assert user_config.get_setting("msdial.exe").value == str(exe)
    use_config(monkeypatch, tmp_path, {"library.lbm": tmp_path / "x.lbm2"})
    assert user_config.get_setting("msdial.exe").value == str(exe)  # 重ねる
    assert user_config.get_setting("library.lbm").value == str(tmp_path / "x.lbm2")
    use_config(monkeypatch, tmp_path, {"msdial.exe": None})
    assert user_config.get_setting("msdial.exe") is None


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
