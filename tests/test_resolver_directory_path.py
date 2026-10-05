"""`file_path` にフォルダが渡されたときのリゾルバの挙動。

Use-LLLM 上のモデルが `load_dataset` と同じフォルダを `curation_review(file_path=...)`
に渡し、フォルダをファイルとして open して `[Errno 13] Permission denied` で落ちた
（2026-10-05）。フォルダはその中の探索先として扱う。
"""

from pathlib import Path

import pytest

from lipidmix.core import mcp_core
from lipidmix.core.path_resolvers import resolve_arf2_file_path, resolve_arf_file_path


@pytest.fixture
def default_dir(tmp_path, monkeypatch):
    directory = tmp_path / "default"
    directory.mkdir()
    (directory / "AlignmentResult_2026_01_01_00_00_00.arf2").write_bytes(b"x")
    monkeypatch.setattr(mcp_core, "DATA_DIR", directory)
    return directory


def test_directory_resolves_to_latest_batch_inside(tmp_path, default_dir):
    target = tmp_path / "neg"
    target.mkdir()
    (target / "AlignmentResult_2026_09_09_17_31_52.arf2").write_bytes(b"x")
    (target / "AlignmentResult_2026_07_21_14_02_08.arf2").write_bytes(b"x")

    resolved = resolve_arf2_file_path(str(target))

    assert Path(resolved) == target / "AlignmentResult_2026_09_09_17_31_52.arf2"


def test_directory_keeps_peakproperties_preference(tmp_path, default_dir):
    target = tmp_path / "neg"
    target.mkdir()
    (target / "AlignmentResult_2026_09_09_17_31_52_DriftSpots.arf").write_bytes(b"x")
    (target / "AlignmentResult_2026_09_09_17_31_52_PeakProperties.arf").write_bytes(b"x")

    resolved = resolve_arf_file_path(str(target))

    assert Path(resolved).name.endswith("_PeakProperties.arf")


def test_directory_without_match_does_not_fall_back_to_default(tmp_path, default_dir):
    # 指定したフォルダに無いのに既定フォルダの別データセットを黙って掴むと、
    # 利用者の意図と違うデータを解析した結果が返る。
    target = tmp_path / "empty"
    target.mkdir()

    assert resolve_arf2_file_path(str(target)) is None


def test_explicit_file_is_still_respected(tmp_path, default_dir):
    target = tmp_path / "picked.arf2"
    target.write_bytes(b"x")

    assert resolve_arf2_file_path(str(target)) == str(target)
