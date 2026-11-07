"""テスト用の設定ファイル（`lipidmix.local.toml`）を tmp に書いて `LIPIDMIX_CONFIG` で指す。

外部資産の場所・置き場所（`[paths]`）・サーバの振る舞い（`[server]`）は設定ファイルだけで
指す（環境変数は読まない。spec 2026-11-06 §7.4）。conftest は `LIPIDMIX_CONFIG` を存在しない
パスへ向けているので、設定が要るテストはここで書く。研究室ライブラリの実パス・実ファイル名は
書かない（tmp の偽物を指す）。

pytest の関数は `use_config(monkeypatch, tmp_path, ...)`、`monkeypatch` の無い
`unittest.TestCase` は `self.enterContext(config_env(self.tmp, ...))`。
"""
from __future__ import annotations

import contextlib
import json
import os
import tomllib
from pathlib import Path
from unittest import mock

from metabolomix.core import user_config

#: 書き出す先のファイル名（`tmp_path` の直下）。
CONFIG_NAME = "lipidmix.local.toml"


def _toml_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    text = str(value)
    if "'" in text or "\n" in text or "\r" in text:
        # リテラル文字列に書けない値だけ基本文字列（JSON の文字列は TOML の基本文字列として読める）。
        return json.dumps(text, ensure_ascii=False)
    return f"'{text}'"  # パスはリテラル文字列（`\` をエスケープとして読ませない）


def _flatten(data: dict) -> dict[str, object]:
    flat: dict[str, object] = {}
    for section, body in data.items():
        if isinstance(body, dict):
            for name, value in body.items():
                flat[f"{section}.{name}"] = value
    return flat


def write_config(directory: Path, values: dict[str, object]) -> Path:
    """`values`（`"<節>.<名前>"` → 値）を `directory/lipidmix.local.toml` に書く（指しはしない）。

    既にあるファイルの値に重ねる（値 None はそのキーを消す）。`Path` は文字列（リテラル文字列
    `'...'`）で書く。相対パスは設定ファイルのフォルダ（`directory`）基準で解かれる。
    """
    path = Path(directory) / CONFIG_NAME
    merged: dict[str, object] = {}
    if path.is_file():
        merged = _flatten(tomllib.loads(path.read_text(encoding="utf-8")))
    for key, value in values.items():
        if value is None:
            merged.pop(key, None)
        else:
            merged[key] = value
    sections: dict[str, list[str]] = {}
    for key, value in merged.items():
        section, name = key.split(".", 1)
        sections.setdefault(section, []).append(f"{name} = {_toml_value(value)}")
    text = "".join(f"[{section}]\n" + "".join(f"{line}\n" for line in lines)
                   for section, lines in sections.items())
    path.parent.mkdir(parents=True, exist_ok=True)
    # 書きかけを読ませない（別スレッドの受付が同時に読むことがある）。置き換えは原子的。
    staging = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    staging.write_text(text, encoding="utf-8")
    os.replace(staging, path)
    return path


def use_config(monkeypatch, tmp_path: Path, values: dict[str, object]) -> Path:
    """`values` を `tmp_path/lipidmix.local.toml` に書き（`write_config`）、そこを指す。

    同じ `tmp_path` で何度呼んでもよい（値を重ねる）。`user_config` の解析キャッシュも消す。
    既定のライブラリのキャッシュ（設定ファイルのフォルダの `data/.library-cache`）も `tmp_path` に入る。
    """
    path = write_config(tmp_path, values)
    monkeypatch.setenv(user_config.CONFIG_ENV, str(path))
    user_config._cache.clear()
    return path


@contextlib.contextmanager
def config_env(directory: Path, values: dict[str, object]):
    """`use_config` の `monkeypatch` を使わない版（`unittest.TestCase` 用）。抜けると元へ戻す。"""
    path = write_config(directory, values)
    with mock.patch.dict(os.environ, {user_config.CONFIG_ENV: str(path)}):
        user_config._cache.clear()
        try:
            yield path
        finally:
            user_config._cache.clear()
