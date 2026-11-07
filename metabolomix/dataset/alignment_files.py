"""解析フォルダから、同じバッチのアラインメント出力一式を 1 組選ぶ。

フォルダ直下だけを見る（サブフォルダは別の解析なので掴まない）。選ばなかった `.arf` は
`not_chosen` に名前で残し、呼び出し側が「他にもあった」と伝えられるようにする。
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.path_resolvers import (
    _pick_latest, _select_latest_batch, _timestamped_first)
from metabolomix.msdial.alignment_spot import alignment_sibling
from metabolomix.msdial.classes import MDDATA_READ_ERRORS, resolve_mddata_path
from metabolomix.msdial.tags import alignment_tag_path

# MS-DIAL が実際に書く綴りは DriftSopts。DriftSpots も受ける。
_ARF_VARIANT_RE = re.compile(r"_(?:PeakProperties|DriftSopts|DriftSpots)\.arf$", re.IGNORECASE)
_PEAK_PROPERTIES_RE = re.compile(r"_PeakProperties\.arf$", re.IGNORECASE)


@dataclass(frozen=True)
class AlignmentSelection:
    arf: Path
    arf2: Path
    dcl: Path | None
    eic: Path | None
    tags: Path | None
    mddata: Path | None
    loaded_libraries: tuple[Path, ...]
    not_chosen: tuple[str, ...]
    #: `.mddata` を決められなかった理由（読み取りの失敗）。無い・決まったときは None。
    mddata_error: str | None = None


def _arf_stem(name: str) -> str:
    match = _ARF_VARIANT_RE.search(name)
    if match:
        return name[: match.start()]
    return name[: -len(".arf")]


def _arf2_stem(name: str) -> str:
    return name[: -len(".arf2")]


def _list_arfs(folder: Path) -> list[Path]:
    return sorted(p for p in folder.iterdir() if p.is_file() and p.name.lower().endswith(".arf"))


def _choose_arf(candidates: list[Path]) -> Path:
    batch = _timestamped_first(_select_latest_batch([str(p) for p in candidates]))
    preferred = [p for p in batch if _PEAK_PROPERTIES_RE.search(p)] or batch
    return Path(_pick_latest(preferred))


def select_alignment(path: str | Path) -> AlignmentSelection:
    path = Path(path)
    lowered = path.name.lower()
    if path.is_dir():
        folder = path
        arfs = _list_arfs(folder)
        if not arfs:
            raise DomainError("ALIGNMENT_NOT_FOUND", f".arf が見つからない: {folder}")
        arf = _choose_arf(arfs)
    elif lowered.endswith(".arf") and path.is_file():
        folder = path.parent
        arfs = _list_arfs(folder)
        arf = path
    elif lowered.endswith(".arf2") and path.is_file():
        folder = path.parent
        arfs = _list_arfs(folder)
        stem = _arf2_stem(path.name).lower()
        same = [p for p in arfs if _arf_stem(p.name).lower() == stem]
        if not same:
            raise DomainError("ALIGNMENT_NOT_FOUND", f"同じバッチの .arf が無い: {path}")
        arf = _choose_arf(same)
    else:
        raise DomainError("ALIGNMENT_NOT_FOUND", f"フォルダか .arf / .arf2 ではない、または存在しない: {path}")

    stem = _arf_stem(arf.name)
    arf2 = folder / f"{stem}.arf2"
    if not arf2.is_file():
        # 大文字小文字違いの拡張子を拾う
        matches = [p for p in folder.iterdir()
                   if p.is_file() and p.name.lower() == f"{stem}.arf2".lower()]
        if not matches:
            raise DomainError(
                "ALIGNMENT_INCOMPLETE", f"同じバッチの .arf2 が無い: {arf.name}",
                {"arf": arf.name},
            )
        arf2 = matches[0]

    tags = alignment_tag_path(arf2)
    mddata_error = None
    try:
        mddata = resolve_mddata_path(arf2)
    except MDDATA_READ_ERRORS as exc:
        # .mddata が読めなくても選択自体は成り立つ。読み込みの注記に残して先へ進む
        # （プログラムの誤りまでは握りつぶさない）。
        mddata = None
        mddata_error = f"{type(exc).__name__}: {exc}"
    loaded = tuple(sorted(folder.glob("*_Loaded.*.dbs")))
    return AlignmentSelection(
        arf=arf,
        arf2=arf2,
        dcl=alignment_sibling(arf2, ".dcl"),
        eic=alignment_sibling(arf2, ".EIC.aef"),
        tags=tags if tags.is_file() else None,
        mddata=mddata,
        loaded_libraries=loaded,
        not_chosen=tuple(p.name for p in arfs if p != arf),
        mddata_error=mddata_error,
    )
