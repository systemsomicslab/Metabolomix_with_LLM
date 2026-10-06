"""雛形 lipidmix.example.toml がそのまま設定ファイルとして読めること。

雛形が壊れていると、外部の利用者は複製した最初の一歩で CONFIG_INVALID に当たる。
"""
from __future__ import annotations

import shutil

from metabolomix.core import user_config


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
