"""library_mode="msp_only" / "msp_then_lbm" のテスト用: tmp に偽の研究室 MSP（と LBM）と設定ファイルを作る。

研究室ライブラリの実パス・実ファイル名は使わない（外部流出禁止の資産）。設定ファイルは
tmp に書いて `LIPIDMIX_CONFIG` で指す（`tests/config_fixtures.py` の `use_config`。
同じ tmp_path に書いた他の設定（`[msdial] exe` など）には値を重ねる）。
"""
from __future__ import annotations

from pathlib import Path

from tests.config_fixtures import use_config


def write_lab_msp_config(tmp_path: Path, monkeypatch, *, positive: bool = True,
                         negative: bool = True, directory_name: str = "lab_library",
                         lbm: bool = False) -> dict[str, Path]:
    """`lbm=True` なら `[library] lbm` に偽の LBM（`lab_lipids.lbm2`）も書く（戻り値の "lbm"）。"""
    lib_dir = tmp_path / directory_name
    lib_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    values: dict[str, object] = {}
    for polarity, enabled, key in (("positive", positive, "msp_positive"),
                                   ("negative", negative, "msp_negative")):
        if not enabled:
            values[f"library.{key}"] = None  # 書き直すたびに [library] を丸ごと置き換える
            continue
        path = lib_dir / f"lab_{polarity[:3]}.msp"
        path.write_text("NAME: fake\nPRECURSORMZ: 100\nNum Peaks: 0\n", encoding="ascii")
        values[f"library.{key}"] = path
        paths[polarity] = path
    values["library.lbm"] = None
    if lbm:
        path = lib_dir / "lab_lipids.lbm2"
        path.write_bytes(b"fake-lab-lbm-library")
        values["library.lbm"] = path
        paths["lbm"] = path
    use_config(monkeypatch, tmp_path, values)
    return paths
