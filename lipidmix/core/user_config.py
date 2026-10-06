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
