"""400 検体規模で最初に当たる失敗（メモリ不足）が理由を持って返ることを固定する。

`str(MemoryError())` は空文字なので、素の `except Exception` に飲まれると
`[ERROR] ARF解析に失敗しました: ` と理由なしで返る。実測モデルは
RSS ≈ 8.0 KB × 行数、行数 = スポット × 注入（`docs/HISTRY.md` 2026-09-20）。
"""
import io
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

import server
from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state


def _group(rows: int = 3) -> list:
    """`extract_arf_data` が受け付ける最小のスポット（行はすべてリスト）。"""
    return [[0] * 24 + ["PC 34:1"] + [0] * 5 for _ in range(rows)]


def test_deserialize_reports_how_far_it_got_when_memory_runs_out(monkeypatch):
    def exploding_groups(data):
        yield 0, 0, _group()
        yield 0, 1, _group()
        raise MemoryError()

    monkeypatch.setattr(arf_reader, "_iter_arf_peak_groups", exploding_groups)

    with pytest.raises(arf_reader.ArfMemoryError) as excinfo:
        arf_reader.deserialize(io.BytesIO(b"ignored"))

    assert isinstance(excinfo.value, MemoryError)
    assert excinfo.value.spots == 2
    assert excinfo.value.rows == 6


class _FakeArfState:
    def __init__(self, error):
        self._error = error
        self.features = None
        self.filtered_features = None
        self.tag_index = {}
        self.class_index = None
        self.excluded_samples = set()
        self.excluded_spots = set()

    def load_data(self, file_path, tag_directory=None):
        raise self._error


class _FakeSession:
    def __init__(self, error):
        self.arf = _FakeArfState(error)

    def maybe_prepend_caveat(self, text, topic=None):
        return text


def _run_arf_parser(error) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "PeakProperties.arf"
        path.touch()
        with patch.object(session_state, "session", _FakeSession(error)):
            return server.arf_parser(str(path))


def test_arf_parser_reports_rows_and_required_ram_on_memory_error():
    result = _run_arf_parser(arf_reader.ArfMemoryError(spots=3125, rows=1_250_000))

    assert "メモリ不足" in result
    assert "3,125 スポット" in result
    assert "1,250,000 行" in result
    assert "10.2 GB" in result


def test_arf_parser_says_the_row_count_is_unknown_when_it_cannot_be_determined():
    result = _run_arf_parser(MemoryError())

    assert "メモリ不足" in result
    assert "不明" in result
    # 理由なしの空メッセージへ戻っていないこと
    assert not result.rstrip().endswith("しました:")
