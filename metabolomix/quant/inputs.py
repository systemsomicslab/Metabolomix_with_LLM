"""定量の入力表（内部標準の参照表・付加体の採用表）の読み込み。

内部標準の参照表は設定ファイルに既定を置かず、呼び出しごとに `reference_path` で指す。
付加体の採用表だけ `[quant] adducts` で指す。戻り値の封筒には、表の sha256・単位・
ファイル名の末尾だけを載せる（パスと中身は載せない）。
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.user_config import get_setting

METHODS = ("cell_coefficient", "standard_concentration")


@dataclass(frozen=True)
class ReferenceInput:
    text: str
    sha256: str
    label: str  # ファイル名の末尾だけ
    path: Path | None
    method: str
    concentrations: dict[str, float] | None
    solvent_volume: float | None
    cells_millions: float | None
    unit: str | None


def quant_error(code: str, message: str, details: dict | None = None) -> str:
    """定量ツールのエラー封筒（`mztab_error` と同じ形）。"""
    return mztab_error(code, message, details)


def _bad(message: str, **details) -> DomainError:
    return DomainError("QUANT_BAD_REQUEST", message, details or None)


def reference_key(label: str) -> str:
    """記録の `inputs` での参照表のキー（絶対パスは載せない。is_stale の現在の入力も同じキーで作る）。"""
    return f"reference:{label}"


def adduct_key(label: str) -> str:
    """記録の `inputs` での付加体の採用表のキー（`reference_key` と同じ流儀）。"""
    return f"adducts:{label}"


def _read(path: Path) -> tuple[str, str]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise DomainError("QUANT_REFERENCE_UNREADABLE",
                          f"表を読めません: {path.name}（{type(exc).__name__}）",
                          {"label": path.name}) from None
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise _bad(f"表を UTF-8 として読めません: {path.name}（UTF-8 で保存し直してください）",
                   label=path.name) from None
    return text, hashlib.sha256(data).hexdigest()


def reference_input(reference_path: str | None, *, method: str = "cell_coefficient",
                    concentrations: dict | None = None,
                    solvent_volume: float | None = None,
                    cells_millions: float | None = None,
                    unit: str | None = None) -> ReferenceInput:
    if not reference_path or not str(reference_path).strip():
        raise _bad("reference_path は必須です（内部標準の参照表の CSV。設定ファイルの既定は読みません）。",
                   field="reference_path")
    if method not in METHODS:
        raise _bad(f"method は {' / '.join(METHODS)} のどちらかです。", field="method")
    if concentrations is not None:
        if method == "standard_concentration":
            raise _bad("濃度方式（standard_concentration）では concentrations を渡せません。"
                       "係数方式の入力です。", field="concentrations")
        if solvent_volume is None or cells_millions is None:
            raise _bad("concentrations を渡すときは solvent_volume と cells_millions も必須です。",
                       field="concentrations")
    elif method == "standard_concentration" and (solvent_volume is not None
                                                 or cells_millions is not None):
        raise _bad("濃度方式（standard_concentration）では solvent_volume・cells_millions を"
                   "渡せません。", field="method")
    path = Path(reference_path).expanduser()
    if not path.is_file():
        raise DomainError("QUANT_REFERENCE_NOT_FOUND",
                          f"参照表が見つかりません: {path.name}", {"label": path.name})
    text, sha = _read(path)
    return ReferenceInput(
        text=text, sha256=sha, label=path.name, path=path, method=method,
        concentrations={str(k): float(v) for k, v in concentrations.items()}
        if concentrations is not None else None,
        solvent_volume=solvent_volume, cells_millions=cells_millions, unit=unit)


def adduct_table() -> ReferenceInput | None:
    """`[quant] adducts` の表。未設定なら None、指す先が無ければ QUANT_REFERENCE_NOT_FOUND。"""
    setting = get_setting("quant.adducts")
    if setting is None:
        return None
    path = Path(setting.value)
    if not path.is_file():
        raise DomainError("QUANT_REFERENCE_NOT_FOUND",
                          f"[quant] adducts の指す先が見つかりません: {path.name}",
                          {"label": path.name, "setting": "quant.adducts"})
    text, sha = _read(path)
    return ReferenceInput(text=text, sha256=sha, label=path.name, path=path,
                          method="cell_coefficient", concentrations=None,
                          solvent_volume=None, cells_millions=None, unit=None)
