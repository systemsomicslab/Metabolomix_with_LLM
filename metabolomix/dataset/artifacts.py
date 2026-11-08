"""`DatasetState` が覚えている付随ファイルの場所を、役割から引く。

探索はしない——読み込みが決めて記録したパス（`artifact_paths` と `source_files`）だけを返す。
役割は旧名と Console の役割名の対（`dcl` / `msms_evidence`、`eic` / `chromatogram`、
`arf2` / `spot_catalog`、`arf` / `peak_matrix_source`）を同じものとして扱う。
stdlib だけ（`core/handles.py` が遅延 import する）。
"""
from __future__ import annotations

from pathlib import Path

__all__ = ["artifact_path", "artifact_paths"]

#: 役割 → 同じものを指す役割名（`artifact_paths` のキー）。
_ALIASES = {
    "dcl": ("dcl", "msms_evidence"),
    "msms_evidence": ("dcl", "msms_evidence"),
    "eic": ("eic", "chromatogram"),
    "chromatogram": ("eic", "chromatogram"),
    "arf2": ("arf2", "spot_catalog"),
    "spot_catalog": ("arf2", "spot_catalog"),
    "arf": ("arf", "peak_matrix_source"),
    "peak_matrix_source": ("arf", "peak_matrix_source"),
}
#: `source_files`（パス → sha256）から拾うときの拡張子。`.arf` は `.arf2` に当たらないよう末尾で見る。
_SUFFIXES = {"arf2": ".arf2", "spot_catalog": ".arf2", "arf": ".arf", "peak_matrix_source": ".arf"}


def artifact_paths(ds, role: str) -> list[Path]:
    """`role` のファイルのうち、今ディスクにあるものを記録の順に返す。"""
    candidates: list[str] = []
    recorded = getattr(ds, "artifact_paths", None) or {}
    for name in _ALIASES.get(role, (role,)):
        candidates.extend(recorded.get(name) or [])
    suffix = _SUFFIXES.get(role)
    if suffix:
        candidates.extend(p for p in (getattr(ds, "source_files", None) or {})
                          if str(p).endswith(suffix))
    found: list[Path] = []
    for text in candidates:
        path = Path(text)
        if path not in found and path.is_file():
            found.append(path)
    return found


def artifact_path(ds, role: str) -> Path | None:
    """`role` のファイル（無ければ None）。複数あれば最初の 1 本。"""
    found = artifact_paths(ds, role)
    return found[0] if found else None
