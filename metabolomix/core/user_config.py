"""外部資産（MS-DIAL Console・脂質ライブラリ・研究室の参照ライブラリ）と置き場所
（レポート・ライブラリのキャッシュ）とサーバの振る舞い（ハンドルの上限・意味論ダイジェスト）の設定。

spec: docs/superpowers/specs/2026-10-06-user-config-file-design.md・2026-11-06 §7.4

値は設定ファイルだけから引く。設定ファイルは `LIPIDMIX_CONFIG` が指すファイル、
無ければ `<repo>/lipidmix.local.toml`（雛形 `lipidmix.example.toml`）。資産の場所・
置き場所を上書きする環境変数は読まない（spec 2026-11-06 §7.4。`LIPIDMIX_CONFIG` は
起動ごとに変わる設定ファイルの場所なので残る）。

呼ばれるたびに読む（パス・更新時刻・大きさが同じなら前回の解析結果を使う）ので、
設定ファイルの書き換えはサーバを再起動しなくても次の呼び出しから効く。

**stdlib だけの leaf**。pipeline の独立 worker も import するので、`mcp_core` /
`session_state` / `metabolomix.tools.*` を import しない（tests/test_user_config.py が縛る）。
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

#: 設定キー（`<節>.<名前>`）→ 値の種類。キーの正準はここ。
#: - `path`: 文字列。設定ファイルのフォルダ基準で解いた絶対パスにする（`~` も展開する）。
#: - `text`: 文字列。前後の空白だけを落とす（URL・トークン）。
#: - `int`: TOML の整数（bool は拒む）。
SETTINGS: dict[str, str] = {
    "msdial.exe": "path",
    "library.msp_positive": "path",
    "library.msp_negative": "path",
    # 脂質ライブラリ（.lbm2）。Console の自動解決（設定の段）と pipeline の
    # library_mode="msp_then_lbm" の両方がここを読む（旧 `[msdial] lbm` はここへ移した）。
    "library.lbm": "path",
    "mslipidmapper.rscript": "path",
    "mslipidmapper.browser": "path",
    "mslipidmapper.mslipidmapper_api_url": "text",
    "mslipidmapper.api_token": "text",
    # 置き場所（置き場所の環境変数から移した。計画 B-2a の決めたこと 2）。
    # reports_dir は解析フォルダ配下 reports/ より先に試す書き込み先。library_cache は
    # 参照ライブラリの SQLite キャッシュ（既定は設定ファイルのフォルダの data/.library-cache）。
    "paths.reports_dir": "path",
    "paths.library_cache": "path",
    # サーバの振る舞い（同じく環境変数から移した）。
    "server.max_dataset_handles": "int",
    "server.caveat_mode": "text",
}

#: 値を選択肢に限る設定 → 選べる値（小文字。大文字小文字は問わず、小文字にそろえて返す）。
#: 他の値は ConfigInvalidError（黙って既定に落とすと、`off` のつもりの打ち間違いが気づかれない）。
CHOICES: dict[str, tuple[str, ...]] = {"server.caveat_mode": ("digest", "off")}

#: 移した旧キー → 新しいキー。設定ファイルに旧キーがあれば ConfigInvalidError
#: （黙って無視すると、`[msdial] lbm` に書いた LBM が Console の解決から消え、
#: 実行体のフォルダの別の LBM へ落ちる）。
MOVED_KEYS: dict[str, str] = {"msdial.lbm": "library.lbm"}

#: パスではない設定（URL・トークン・整数）。設定ファイル基準のパスに解かない。
NON_PATH_SETTINGS: frozenset[str] = frozenset(
    key for key, kind in SETTINGS.items() if kind != "path")

#: 極性 → 研究室参照ライブラリ（.msp）の設定キー。正準はここ
#: （`core/path_resolvers.LIBRARY_SETTING_KEYS` と `console/method_file.resolve_msp` が参照する。
#: pipeline の worker も使うので、mcp_core を引き込まない leaf のここに置く）。
MSP_SETTING_KEYS: dict[str, str] = {
    "positive": "library.msp_positive",
    "negative": "library.msp_negative",
}

#: 脂質ライブラリ（.lbm2）の設定キー（`[library] lbm`）。
LIBRARY_LBM_SETTING_KEY = "library.lbm"

# 二重引用符の中の `\n` `\t` `\b` `\f` `\r` は TOML の正しいエスケープなので構文エラーに
# ならない。`"C:\new\tool.exe"` は黙って改行とタブ入りの値になる。値の側で拾う。
_CONTROL_CHARS = frozenset("\n\t\b\f\r")
_QUOTE_HINT = ("Windows のパスは単一引用符で囲んでください"
               "（例: exe = 'C:\\MS-DIAL\\MSDIALCUI.exe'）。")


@dataclass(frozen=True)
class Setting:
    """1 つの設定値と、それを書いた設定ファイル。"""

    key: str
    #: 種類 `path` は絶対パス、`text` は前後の空白を除いた文字列、`int` は整数。
    value: str | int
    config_file: str


class ConfigInvalidError(Exception):
    """設定ファイルが読めない（構文・型・壊れたパス・移した旧キー）。メッセージは利用者向け。"""

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
    values: dict[str, str | int]
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


def _section_and_name(key: str) -> tuple[str, str]:
    section, name = key.split(".", 1)
    return section, name


def config_key(key: str) -> str:
    """設定ファイルでの書き方（`[節] 名前`）。"""
    section, name = _section_and_name(key)
    return f"[{section}] {name}"


def _value(path: Path, key: str, value: object) -> str | int | None:
    """設定ファイルの 1 つの値を種類に合わせて検査・整形する。空の文字列は None（未設定）。"""
    kind = SETTINGS[key]
    where = f"{path.name} の {config_key(key)}"
    if kind == "int":
        if type(value) is not int:  # bool は int の部分型なので isinstance では拒めない
            raise ConfigInvalidError(
                f"{where} は整数で書いてください（受け取った型: {type(value).__name__}）。",
                config_file=str(path))
        return value
    if not isinstance(value, str):
        what = "文字列（パス）" if kind == "path" else "文字列"
        raise ConfigInvalidError(
            f"{where} は{what}で書いてください（受け取った型: {type(value).__name__}）。"
            f"{_QUOTE_HINT if kind == 'path' else ''}",
            config_file=str(path))
    if _CONTROL_CHARS & set(value):
        raise ConfigInvalidError(
            f"{where} に改行やタブが入っています。二重引用符の中の \\n や \\t が"
            f"エスケープとして読まれた可能性があります。{_QUOTE_HINT}",
            config_file=str(path))
    value = value.strip()
    if not value:
        return None
    if kind == "path":
        return _resolve(value, path.parent)
    choices = CHOICES.get(key)
    if choices is not None:
        if value.lower() not in choices:
            raise ConfigInvalidError(
                f"{where} は {' / '.join(choices)} のどれかで書いてください（受け取った値: {value!r}）。",
                config_file=str(path))
        return value.lower()
    return value


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

    values: dict[str, str | int] = {}
    unknown: list[str] = []
    for section, body in data.items():
        if not isinstance(body, dict):
            unknown.append(section)
            continue
        for name, raw in body.items():
            key = f"{section}.{name}"
            if key in MOVED_KEYS:
                raise ConfigInvalidError(
                    f"{path.name} の {config_key(key)} は {config_key(MOVED_KEYS[key])} へ移しました。"
                    f"{config_key(MOVED_KEYS[key])} に書き直してください（MS-DIAL Console の"
                    "脂質ライブラリの自動解決も pipeline の library_mode='msp_then_lbm' も"
                    f" {config_key(MOVED_KEYS[key])} だけを読みます）。",
                    config_file=str(path))
            if key not in SETTINGS:
                unknown.append(key)
                continue
            value = _value(path, key, raw)
            if value is not None:
                values[key] = value

    parsed = _Parsed(values=values, unknown_keys=tuple(sorted(unknown)))
    _cache.clear()
    _cache[cache_key] = parsed
    return parsed


def get_setting(key: str) -> Setting | None:
    """`key` の値を設定ファイルから引く。書かれていなければ（設定ファイルが無ければ）None。

    設定ファイルが読めなければ ConfigInvalidError。
    """
    if key not in SETTINGS:
        raise KeyError(key)
    path = config_file_path()
    parsed = _parse(path)
    if parsed is None or key not in parsed.values:
        return None
    return Setting(key=key, value=parsed.values[key], config_file=str(path))


def setting_label(setting: Setting) -> str:
    """エラー文で値の出どころを言う句（「lipidmix.local.toml の [msdial] exe」）。"""
    return f"{Path(setting.config_file).name} の {config_key(setting.key)}"


def _unknown_keys_quietly(path: Path) -> tuple[tuple[str, ...], str | None]:
    """(未知のキー, 読めなかった理由)。エラー文の材料なので送出しない。"""
    try:
        parsed = _parse(path)
    except ConfigInvalidError as exc:
        return (), exc.message
    return (parsed.unknown_keys if parsed else ()), None


def _instruction(key: str, what: str) -> str:
    """「〜を書いてください」を値の種類に合わせる（パスでない設定に「パス」と言わない）。"""
    kind = SETTINGS[key]
    if kind == "path":
        return f"{what}のパスを書いてください"
    if kind == "int":
        return f"{what}を整数で書いてください"
    choices = CHOICES.get(key)
    if choices is not None:
        return f"{what}を {' / '.join(choices)} のどれかで書いてください"
    return f"{what}を書いてください"


def missing_hint(key: str, what: str) -> str:
    """未設定のとき利用者に伝える文。`what` は「MS-DIAL Console の実行体」など。"""
    path = config_file_path()
    instruction = _instruction(key, what)
    if path.is_file():
        head = f"{path.name} の {config_key(key)} に{instruction}"
    else:
        head = (f"{EXAMPLE_FILENAME} を {path} として複製し、"
                f"{config_key(key)} に{instruction}")
    quote = "Windows のパスは単一引用符で囲む。" if SETTINGS[key] == "path" else ""
    message = f"{head}（{quote}書き換えは次の呼び出しから効き、サーバの再起動は要りません）。"
    unknown, _ = _unknown_keys_quietly(path)
    if unknown:
        message += f" 認識できないキーがあります: {', '.join(unknown)}。"
    return message


def describe_missing(key: str, setting: Setting | None = None) -> dict:
    """未設定・指す先が無いエラーの `details` に足す内容。値そのものは載せない。"""
    path = config_file_path()
    details: dict = {
        "setting": key,
        "config_key": config_key(key),
        "source": "config_file" if setting is not None else None,
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
