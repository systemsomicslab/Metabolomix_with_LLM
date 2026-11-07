"""外部資産の場所の環境変数（spec 2026-11-06 §7.4 の 9 つ）を `metabolomix/` が読まない・案内しない。

場所は `lipidmix.local.toml` だけで指す。エラーコード（`MSDIAL_EXE_NOT_FOUND` など）は
名前が違うので当たらない。
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERN = re.compile(
    r"\b(MSDIAL_EXE|MSDIAL_LBM|MSDIAL_MSP_POS|MSDIAL_MSP_NEG|MSDIAL_LIBRARY_LBM"
    r"|MSLIPIDMAPPER_[A-Z_]+)\b(?!_NOT)")


def _hits(text: str) -> list[str]:
    return [m.group(0) for m in PATTERN.finditer(text)]


def test_pattern_skips_error_codes():
    assert _hits("MSDIAL_EXE_NOT_FOUND MSDIAL_EXE_NOT_CONSOLE MSDIAL_EXE_PROBE_TIMEOUT "
                 "MSDIAL_EXECUTION_FAILED") == []
    assert _hits("os.environ['MSDIAL_EXE'] getenv('MSLIPIDMAPPER_API_TOKEN') "
                 "環境変数 MSDIAL_LBM") == ["MSDIAL_EXE", "MSLIPIDMAPPER_API_TOKEN", "MSDIAL_LBM"]


def test_no_asset_env_vars():
    found: list[str] = []
    for path in sorted((ROOT / "metabolomix").rglob("*.py")):
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for name in _hits(line):
                found.append(f"{path.relative_to(ROOT).as_posix()}:{lineno}: {name}")
    assert found == []
