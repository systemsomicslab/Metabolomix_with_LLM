"""400 検体規模で最初に当たる失敗（メモリ不足）が理由を持って返ることを固定する。

`str(MemoryError())` は空文字なので、素の `except Exception` に飲まれると
理由なしで返る（旧 `arf_parser` の失敗文言）。実測モデルは
RSS ≈ 8.0 KB × 行数、行数 = スポット × 注入（`docs/HISTRY.md` 2026-09-20）。
"""
import io

import pytest

from metabolomix.arf import reader as arf_reader


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
