"""描画のロック: ワーカースレッドの matplotlib 描画は PLOT_LOCK 1 本で直列にする。"""
from __future__ import annotations

import re
import threading
import time
from pathlib import Path

import matplotlib.pyplot as plt

from metabolomix.plots import mirror, render
from metabolomix.plots import volcano as volcano_plot

ROOT = Path(__file__).resolve().parent.parent / "metabolomix"
PNG_MAGIC = b"\x89PNG"


def _volcano_payload():
    return volcano_plot.build_volcano_plot_payload({
        "kind": "two_group", "a": "A", "b": "B", "n_a": 3, "n_b": 3,
        "q_threshold": 0.05, "log2fc_threshold": 1.0,
        "volcano": [
            {"feature": "PC 34:1", "log2fc": 1.8, "neg_log10_p": 3.4, "sig": "up"},
            {"feature": "PE 36:2", "log2fc": -2.1, "neg_log10_p": 4.0, "sig": "down"},
            {"feature": "TG 52:3", "log2fc": 0.2, "neg_log10_p": 0.5, "sig": "ns"},
        ]})


def _mirror_payload():
    return mirror.build_mirror_payload(
        [[100.0, 999.0], [200.0, 500.0]], [[100.0, 999.0]],
        [{"mz": 100.0, "measured": 1.0, "reference": 1.0, "matched": True}], title="t")


def test_concurrent_renders_are_serialized():
    events: list[tuple[str, int]] = []

    def traced(fn):
        def wrapper(*args, **kwargs):
            events.append(("in", threading.get_ident()))
            try:
                time.sleep(0.005)  # 重なれば入退出の列が交互にならない
                return fn(*args, **kwargs)
            finally:
                events.append(("out", threading.get_ident()))
        return wrapper

    volcano_build = traced(volcano_plot.render_volcano_plot)
    mirror_build = traced(mirror._build_mirror_figure)
    v_payload, m_payload = _volcano_payload(), _mirror_payload()
    results, errors = [], []

    def work(i):
        try:
            if i % 2:
                results.append(render.render_png(mirror_build, m_payload))
            else:
                results.append(render.render_png(volcano_build, v_payload))
        except BaseException as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=work, args=(i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors
    assert len(results) == 8 and all(r.startswith(PNG_MAGIC) for r in results)
    assert plt.get_fignums() == []
    assert len(events) == 16
    for i in range(0, 16, 2):
        assert events[i][0] == "in" and events[i + 1][0] == "out"
        assert events[i][1] == events[i + 1][1]


def test_render_mirror_keeps_its_bytes_contract():
    assert mirror.render_mirror(_mirror_payload()).startswith(PNG_MAGIC)
    assert plt.get_fignums() == []


def test_render_files_writes_each_format(tmp_path):
    paths = {"png": tmp_path / "a.png", "svg": tmp_path / "a.svg"}
    render.render_files(volcano_plot.render_volcano_plot, paths, _volcano_payload())
    assert paths["png"].read_bytes().startswith(PNG_MAGIC)
    assert b"<svg" in paths["svg"].read_bytes()
    assert plt.get_fignums() == []


def test_render_files_atomic_leaves_no_tmp_and_keeps_old_file_on_failure(tmp_path):
    target = tmp_path / "a.png"
    target.write_bytes(b"old")

    def boom():
        raise RuntimeError("x")

    try:
        render.render_files(boom, {"png": target}, atomic=True)
    except RuntimeError:
        pass
    assert target.read_bytes() == b"old"
    render.render_files(volcano_plot.render_volcano_plot, {"png": target},
                        _volcano_payload(), atomic=True)
    assert target.read_bytes().startswith(PNG_MAGIC)
    assert [p.name for p in tmp_path.iterdir()] == ["a.png"]


# CLI 用の `plot_*`（python -m metabolomix.arf.reader）は MCP のワーカー経路ではない。
_EXCLUDED = {Path("arf") / "reader.py"}


def test_no_direct_figure_to_png_outside_render():
    pattern = re.compile(r"\bfigure_to_png\(|\.savefig\(")
    offenders = []
    for path in ROOT.rglob("*.py"):
        rel = path.relative_to(ROOT)
        if rel == Path("plots") / "render.py" or rel in _EXCLUDED:
            continue
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if pattern.search(line):
                offenders.append(f"{rel}:{n}")
    assert offenders == []
