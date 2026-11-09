"""細かい不具合（spec 2026-10-22 §4 の 9 (a)(b)）。

(a) Bruker の `.d` の中に baf2sql が `analysis.sqlite` を書いても、入力の指紋が変わらない
    （Console 自身の読み込みで INPUT_CHANGED にならない）。
(b) メソッドファイルに ASCII 以外があれば、計画の前に行と文字を示して止める
    （`'ascii' codec` の例外にしない）。
"""
from __future__ import annotations

import json
import os
import time

import pytest

from metabolomix.console import execution, job_manager
from metabolomix.console import method_file as method_file_mod
from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import inputs
from tests.config_fixtures import use_config


def _bruker(root):
    d = root / "S1.d"
    d.mkdir(parents=True)
    (d / "analysis.baf").write_bytes(b"baf")
    return d


def _touch_cache(d):
    time.sleep(0.02)
    (d / "analysis.sqlite").write_bytes(b"cache")


def test_vendor_cache_is_recognised_only_inside_d(tmp_path):
    assert job_manager.is_vendor_cache(tmp_path / "S1.d" / "analysis.sqlite")
    assert job_manager.is_vendor_cache(tmp_path / "S1.D" / "Analysis.sqlite-journal")
    assert not job_manager.is_vendor_cache(tmp_path / "S1.raw" / "analysis.sqlite")
    assert not job_manager.is_vendor_cache(tmp_path / "S1.d" / "analysis.baf")


def test_execution_fingerprint_ignores_cache_written_into_d(tmp_path):
    d = _bruker(tmp_path)
    before = execution.raw_fingerprint([str(d)])
    _touch_cache(d)
    assert execution.raw_fingerprint([str(d)]) == before
    time.sleep(0.02)
    (d / "analysis.baf").write_bytes(b"changed")
    assert execution.raw_fingerprint([str(d)]) != before


def test_execution_fingerprint_unchanged_for_plain_files(tmp_path):
    raw = tmp_path / "a.mzML"
    raw.write_bytes(b"x")
    st = os.stat(raw)
    from metabolomix.core.atomic_io import canonical_hash
    assert execution.raw_fingerprint([str(raw)]) == canonical_hash(
        [{"path": str(raw), "size": st.st_size, "mtime_ns": st.st_mtime_ns}])


def test_verify_inputs_ignores_cache_in_old_snapshots(tmp_path):
    d = _bruker(tmp_path)
    _touch_cache(d)
    stats = [{"relative_path": rel, "size": (tmp_path / rel).stat().st_size,
              "mtime_ns": (tmp_path / rel).stat().st_mtime_ns, "role": "primary"}
             for rel in ("S1.d/analysis.baf", "S1.d/analysis.sqlite")]
    snapshot = {"source_root": str(tmp_path), "raw_stat": stats, "staged_files": []}
    time.sleep(0.02)
    (d / "analysis.sqlite").write_bytes(b"cache rewritten by baf2sql")
    inputs.verify_inputs(snapshot)                     # INPUT_CHANGED にならない
    time.sleep(0.02)
    (d / "analysis.baf").write_bytes(b"really changed")
    with pytest.raises(DomainError) as info:
        inputs.verify_inputs(snapshot)
    assert info.value.code == "INPUT_CHANGED"


def _method(tmp_path, text):
    path = tmp_path / "method.txt"
    path.write_bytes(text.encode("utf-8"))
    return path


def test_check_ascii_names_line_and_character(tmp_path):
    path = _method(tmp_path, "Ion mode: Positive\nMemo: 日本語のメモ書き\nMS1 tolerance: 0.01\n")
    with pytest.raises(DomainError) as info:
        method_file_mod.check_ascii(path)
    assert info.value.code == "METHOD_FILE_NOT_ASCII"
    assert (info.value.details["line"], info.value.details["column"]) == (2, 7)
    assert info.value.details["char"] == "日" and info.value.details["count"] == 8
    method_file_mod.check_ascii(_method(tmp_path, "Ion mode: Positive\n"))


def test_write_effective_method_file_refuses_non_ascii_without_crashing(tmp_path):
    src = _method(tmp_path, "Ion mode: Positive\nMemo: é\n")
    with pytest.raises(DomainError) as info:
        method_file_mod.write_effective_method_file(src, tmp_path / "out.txt",
                                                    {"Ion mode": "Negative"})
    assert info.value.code == "METHOD_FILE_NOT_ASCII"


def test_non_ascii_in_an_overridden_line_is_fine(tmp_path):
    """上書きで置き換わる行の ASCII 以外は効かない（以前も通っていた。挙動を変えない）。"""
    src = _method(tmp_path, "Ion mode: ポジティブ\nMS1 tolerance: 0.01\n")
    out = method_file_mod.write_effective_method_file(src, tmp_path / "out.txt",
                                                      {"Ion mode": "Negative"})
    assert out.read_text(encoding="ascii").splitlines()[0] == "Ion mode: Negative"
