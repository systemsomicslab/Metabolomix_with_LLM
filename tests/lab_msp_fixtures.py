"""library_mode="msp_only" / "msp_then_lbm" のテスト用: tmp に偽の研究室 MSP（と LBM）と設定ファイルを作る。

研究室ライブラリの実パス・実ファイル名は使わない（外部流出禁止の資産）。設定ファイルは
tmp に書いて `LIPIDMIX_CONFIG` で指す（tests/conftest.py の隔離 fixture の流儀）。
"""
from __future__ import annotations

from pathlib import Path


def write_lab_msp_config(tmp_path: Path, monkeypatch, *, positive: bool = True,
                         negative: bool = True, directory_name: str = "lab_library",
                         lbm: bool = False) -> dict[str, Path]:
    """`lbm=True` なら `[library] lbm` に偽の LBM（`lab_lipids.lbm2`）も書く（戻り値の "lbm"）。"""
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
    if lbm:
        path = lib_dir / "lab_lipids.lbm2"
        path.write_bytes(b"fake-lab-lbm-library")
        lines.append(f"lbm = '{path}'")
        paths["lbm"] = path
    config =tmp_path / "lipidmix.local.toml"
    config.write_text("\n".join(lines) + "\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(config))
    return paths


def remember_library_for_mztab(monkeypatch, msp_path: Path) -> None:
    """mzTab-M を直接読んだデータセットが、`msp_path` を同定に使ったと覚えているようにする。

    mzTab-M からは同定に使ったライブラリが分からない（`unknown`）ので、`dataset_handle_open
    (library="auto")` が記憶から開く経路を、解析フォルダを組まずに通すための差し込み。
    """
    from metabolomix.dataset import library_memory
    from metabolomix.handoff.schema import sha256_file

    msp_path = Path(msp_path)
    record = {"status": "known", "kind": "msp", "file": msp_path.name,
              "path": str(msp_path.resolve()), "sha256": sha256_file(msp_path),
              "origin": "effective_method", "candidates": [], "lbm_also": False}
    monkeypatch.setattr(library_memory, "unknown_record", lambda: dict(record))
