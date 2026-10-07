"""library_mode="msp_only" のテスト用: tmp に偽の研究室 MSP と設定ファイルを作る。

研究室ライブラリの実パス・実ファイル名は使わない（外部流出禁止の資産）。設定ファイルは
tmp に書いて `LIPIDMIX_CONFIG` で指す（tests/conftest.py の隔離 fixture の流儀）。
"""
from __future__ import annotations

from pathlib import Path


def write_lab_msp_config(tmp_path: Path, monkeypatch, *, positive: bool = True,
                         negative: bool = True, directory_name: str = "lab_library") -> dict[str, Path]:
    lib_dir = tmp_path / directory_name
    lib_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    lines = ["[library]"]
    for polarity, enabled, key in (("positive", positive, "msp_positive"),
                                   ("negative", negative, "msp_negative")):
        if not enabled:
            continue
        path = lib_dir / f"lab_{polarity[:3]}.msp"
        path.write_text("NAME: fake\nPRECURSORMZ: 100\nNum Peaks: 0\n", encoding="ascii")
        lines.append(f"{key} = '{path}'")
        paths[polarity] = path
    config = tmp_path / "lipidmix.local.toml"
    config.write_text("\n".join(lines) + "\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(config))
    return paths
