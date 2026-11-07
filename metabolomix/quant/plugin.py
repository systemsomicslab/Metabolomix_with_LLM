"""任意の規則プラグイン（パッケージ `msdial_rules`）の読み込み。

ms-data-parser は規則を依存として宣言しない。`lipidmix.local.toml` の `[quant] rules` が
パッケージのフォルダを指すときだけ、固有の名前で読み込む（`sys.modules` に `msdial_rules`
を残さない）。面の仕様は docs/quant-rules-plugin.md。
"""
from __future__ import annotations

import hashlib
import importlib
import importlib.util
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.user_config import get_setting, missing_hint

SUPPORTED_API_VERSION = 1
UNAVAILABLE = "quant_rules_unavailable"
INCOMPATIBLE = "quant_rules_incompatible"
REJECTED = "QUANT_RULES_REJECTED"


@dataclass(frozen=True)
class RulesPlugin:
    api: ModuleType
    rules_version: str
    digest: str
    label: str


#: (解決済みフォルダ, digest) → プラグイン
_cache: dict[tuple[str, str], RulesPlugin] = {}


def _digest(folder: Path) -> str:
    h = hashlib.sha256()
    files = [p for p in folder.rglob("*.py") if "__pycache__" not in p.parts]
    for path in sorted(files, key=lambda p: p.relative_to(folder).as_posix()):
        h.update(path.relative_to(folder).as_posix().encode("utf-8") + b"\0")
        h.update(path.read_bytes() + b"\0")
    return h.hexdigest()


_MAX_TEXT = 300
_STOP = r"""[^\s'"<>|*?:,;)\]]"""
_PATH = re.compile(r"(?:[A-Za-z]:[\\/]|\\\\)" + _STOP + r"*|(?<![\w.:/])/" + _STOP + "+")


def _tail(match: re.Match) -> str:
    parts = [p for p in re.split(r"[\\/]+", match.group(0)) if p]
    return parts[-1] if parts else ""


def sanitize_text(text) -> str:
    """プラグイン由来の文字列を封筒へ写す前に整える: 絶対パスはファイル名の末尾へ、300 字まで。"""
    text = _PATH.sub(_tail, str(text))
    return text if len(text) <= _MAX_TEXT else text[:_MAX_TEXT] + "..."


def sanitize_details(details) -> dict:
    """JSON のスカラ値（str は整える・int・float・bool・None）の値だけ残す。"""
    out = {}
    if not isinstance(details, dict):
        return out
    for key, value in details.items():
        if isinstance(value, str):
            out[str(key)] = sanitize_text(value)
        elif value is None or isinstance(value, (bool, int, float)):
            out[str(key)] = value
    return out


def _unavailable(message: str) -> DomainError:
    return DomainError(UNAVAILABLE, sanitize_text(message))


def _import(folder: Path, digest: str) -> ModuleType:
    name = f"_lipidmix_quant_rules_{digest[:8]}"
    spec = importlib.util.spec_from_file_location(
        name, folder / "__init__.py", submodule_search_locations=[str(folder)])
    if spec is None or spec.loader is None:
        raise _unavailable("[quant] rules のフォルダを Python のパッケージとして読めません。")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
        return importlib.import_module(f"{name}.api")
    except BaseException:
        for key in [k for k in sys.modules if k == name or k.startswith(name + ".")]:
            del sys.modules[key]
        raise


def load_rules() -> RulesPlugin:
    setting = get_setting("quant.rules")
    if setting is None:
        raise _unavailable(
            "規則プラグインが未設定です。"
            + missing_hint("quant.rules", "規則パッケージ（msdial_rules）のフォルダ"))
    folder = Path(setting.value)
    if not folder.is_dir():
        raise _unavailable(f"[quant] rules の指すフォルダが見つかりません: {folder.name}")
    if not (folder / "__init__.py").is_file():
        raise _unavailable(f"[quant] rules のフォルダに __init__.py がありません: {folder.name}")
    digest = _digest(folder)
    key = (str(folder.resolve()), digest)
    cached = _cache.get(key)
    if cached is not None:
        return cached
    try:
        api = _import(folder, digest)
    except Exception as exc:
        raise _unavailable(
            f"規則パッケージを import できません: {type(exc).__name__}: {exc}") from exc
    version = getattr(api, "API_VERSION", None)
    if version != SUPPORTED_API_VERSION:
        raise DomainError(
            INCOMPATIBLE,
            f"規則プラグインの API_VERSION は {version!r} ですが、この ms-data-parser は "
            f"{SUPPORTED_API_VERSION} だけを扱います。",
            {"api_version": version, "supported": SUPPORTED_API_VERSION})
    plugin = RulesPlugin(api=api, rules_version=str(getattr(api, "RULES_VERSION", "")),
                         digest=digest, label=folder.name)
    for stale in [k for k in _cache if k[0] == key[0]]:
        del _cache[stale]
    _cache[key] = plugin
    return plugin


def rules_settings(plugin: RulesPlugin) -> dict:
    """記録の `settings.rules` に入れる。"""
    return {"rules_version": plugin.rules_version, "digest": plugin.digest,
            "api_version": plugin.api.API_VERSION}


def call(plugin: RulesPlugin, name: str, **kwargs):
    try:
        return getattr(plugin.api, name)(**kwargs)
    except ValueError as exc:
        code = getattr(exc, "code", None)
        details = sanitize_details(getattr(exc, "details", None))
        raise DomainError(REJECTED, sanitize_text(exc),
                          {"rule_code": sanitize_text(code) if isinstance(code, str) else None,
                           **details}) from exc
