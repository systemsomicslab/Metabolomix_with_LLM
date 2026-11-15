"""raw の取得時刻を読む道具をイベントループから外す（レビュー指摘 R17）。

`list_data_files(all_files=True)` と `pipeline_sample_sheet` は raw ごとに見出しを開く（.wiff の OLE・Thermo の
見出し・Bruker の sqlite・mzML）。数百のファイルやネットワークの共有では数十秒かかりうるので、
`handle_dispatch.ALWAYS_THREADED` に載せてワーカースレッドで回す（その間も `job_status` などが答えられる）。
見出しの読み取りを遅くして、呼び出しの間にイベントループが進むかを数える。fixture はテスト自身が tmp に作る。
"""
from __future__ import annotations

import asyncio
import json
import time
from datetime import timedelta

import pytest

import server
from metabolomix.core import handle_dispatch
from metabolomix.rawmeta import acquisition
from tests.test_rawmeta import T0, _isa, write_mzml

_DELAY = 0.1   # raw 1 本あたりに足す秒数


@pytest.fixture
def slow_reads(monkeypatch):
    real = acquisition.read_acquisition_times

    def slow(path):
        time.sleep(_DELAY)
        return real(path)

    monkeypatch.setattr(acquisition, "read_acquisition_times", slow)


def _ticks_during(name: str, args: dict) -> tuple[int, str]:
    """`name` を呼ぶ間に 20 ms ごとの ticker が何回進んだかと、戻り値の本文。"""
    ticks = 0

    async def main():
        nonlocal ticks
        done = asyncio.Event()

        async def ticker():
            nonlocal ticks
            while not done.is_set():
                await asyncio.sleep(0.02)
                ticks += 1

        task = asyncio.create_task(ticker())
        try:
            return await server.mcp.call_tool(name, args)
        finally:
            done.set()
            await task

    content = asyncio.run(main())
    blocks = content[0] if isinstance(content, tuple) else content
    return ticks, blocks[0].text


@pytest.mark.parametrize("name", ["list_data_files", "pipeline_sample_sheet"])
def test_raw_scanning_tools_are_threaded(name):
    assert name in handle_dispatch.ALWAYS_THREADED
    assert server.mcp._tool_manager._tools[name].is_async


def test_list_data_files_does_not_block_the_loop(tmp_path, slow_reads):
    for i in range(6):
        write_mzml(tmp_path / f"r{i:02d}.mzML", T0 + timedelta(minutes=10 * i))
    ticks, text = _ticks_during("list_data_files", {"directory": str(tmp_path), "all_files": True})
    # 6 本 × 0.1 秒の間、ループが止まっていれば 1 回ほどしか進まない。
    assert ticks >= 10, f"一覧の間にイベントループが進んでいない（{ticks} 回）"
    assert "r05.mzML" in text


def test_pipeline_sample_sheet_does_not_block_the_loop(tmp_path, slow_reads):
    raw, s_table, a_table = _isa(tmp_path)
    ticks, text = _ticks_during("pipeline_sample_sheet", {
        "raw_dir": str(raw), "sample_table": str(s_table), "factor_columns": ["Genotype"],
        "assay_table": str(a_table)})
    assert ticks >= 10, f"シートを作る間にイベントループが進んでいない（{ticks} 回）"
    assert json.loads(text)["status"] == "success"
