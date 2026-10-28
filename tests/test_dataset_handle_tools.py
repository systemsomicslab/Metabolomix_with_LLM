"""既存のツールの `dataset_handle`（spec 2026-10-22 §4 の 1）。

- 閉じた・無いハンドルは封筒で答え、対話のセッションへ落ちない。
- ハンドルの上での操作は対話のセッション（と別のハンドル）を変えない。
- ハンドルを渡さない呼び出しの戻り値は変更前と同じ（変更前のコードで採った値と比べる）。
"""
from __future__ import annotations

import asyncio
import json
import threading
import time
from pathlib import Path

import pytest

import server
from metabolomix.core import handles, session_state
from metabolomix.library import tools as library_tools
from metabolomix.mztab import identity as mztab_identity
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools, mztab_tools
from metabolomix.tools import reanalysis_tools as R
from tests.reanalysis_fixtures import NEG_FEATURES, write_mztab

GOLDEN = Path(__file__).parent / "fixtures" / "handle_golden_no_handle.json"
HC, RA = ["HC1", "HC2", "HC3"], ["RA1", "RA2", "RA3"]
UNKNOWN = "dh_000000000000"


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


@pytest.fixture
def pos_mztab(tmp_path):
    return write_mztab(tmp_path / "AlignResult-1.mzTab")


def _open(path) -> str:
    out = json.loads(handle_tools.dataset_handle_open(mztab_path=str(path), library="none"))
    return out["dataset_handle"]


def _tools(tmp_path):
    """ハンドルを受けるツール（既存の 13 と dataset_alignment_qc）と最小の引数。"""
    return [
        (mztab_tools.dataset_status, {}),
        (T.dataset_preprocess, {}),
        (T.dataset_pca, {}),
        (T.dataset_differential, {"group_a": HC, "group_b": RA}),
        (T.dataset_export_differential, {"output_path": str(tmp_path / "x.tsv")}),
        (T.dataset_export_evidence_bundle, {"output_dir": str(tmp_path / "b"), "result_id": "r"}),
        (T.dataset_set_sample_metadata, {"manifest_path": str(tmp_path / "m.tsv")}),
        (T.dataset_statistic, {"specification": {}, "matrix_result_id": "m"}),
        (T.dataset_build_matrix, {"recipe": {}}),
        (T.dataset_feature_search, {"feature_ids": ["1"]}),
        (library_tools.library_load, {"ion_mode": "positive"}),
        (library_tools.library_match_feature, {"precursor_mz": 104.107}),
        (library_tools.library_plot_mirror, {}),
        (R.dataset_alignment_qc, {}),
    ]


def test_every_listed_tool_accepts_dataset_handle(tmp_path):
    import inspect
    for tool, _ in _tools(tmp_path):
        assert "dataset_handle" in inspect.signature(tool).parameters, tool.__name__
        assert inspect.signature(tool).parameters["dataset_handle"].default is None


def _fake_pipeline_run(tmp_path, *names: str) -> Path:
    """`input/` に `.arf2`・`.dcl` を置いた pipeline の run の殻（中身は読まれない）。"""
    run = tmp_path / "pipeline_x"
    (run / "input").mkdir(parents=True)
    for name in names:
        (run / "input" / name).write_bytes(b"")
    return run


def test_handle_open_finds_pipeline_arf2(tmp_path, monkeypatch, pos_mztab):
    """pipeline の run では mzTab は console/ の下、`.arf2`・`.dcl` は input/ にある。"""
    from metabolomix.mztab.loading import load_dataset_state

    run = _fake_pipeline_run(tmp_path, "AlignResult-1.arf2", "AlignResult-1.dcl")
    monkeypatch.setattr(handle_tools, "_load",
                        lambda kind, path, allow: load_dataset_state(mztab_path=pos_mztab))

    out = json.loads(handle_tools.dataset_handle_open(pipeline_path=str(run), library="none"))

    assert out["status"] == "success", out
    assert out["arf2"] == "AlignResult-1.arf2"
    assert out["dcl"] == "AlignResult-1.dcl"


def test_handle_open_picks_newest_arf2_name_in_pipeline_input(tmp_path, monkeypatch, pos_mztab):
    from metabolomix.mztab.loading import load_dataset_state

    run = _fake_pipeline_run(tmp_path, "AlignResult-1.arf2", "AlignResult-2.arf2",
                             "AlignResult-2.dcl")
    monkeypatch.setattr(handle_tools, "_load",
                        lambda kind, path, allow: load_dataset_state(mztab_path=pos_mztab))

    out = json.loads(handle_tools.dataset_handle_open(pipeline_path=str(run), library="none"))

    assert out["arf2"] == "AlignResult-2.arf2"
    assert out["dcl"] == "AlignResult-2.dcl"


def test_handle_open_without_pipeline_arf2_stays_null(tmp_path, monkeypatch, pos_mztab):
    from metabolomix.mztab.loading import load_dataset_state

    run = _fake_pipeline_run(tmp_path)
    monkeypatch.setattr(handle_tools, "_load",
                        lambda kind, path, allow: load_dataset_state(mztab_path=pos_mztab))

    out = json.loads(handle_tools.dataset_handle_open(pipeline_path=str(run), library="none"))

    assert out["arf2"] is None and out["dcl"] is None


def test_unknown_handle_is_an_error_not_the_global_session(tmp_path, pos_mztab):
    mztab_tools.dataset_load(mztab_path=str(pos_mztab))       # 対話のセッションにはデータがある
    for tool, kwargs in _tools(tmp_path):
        out = json.loads(tool(**kwargs, dataset_handle=UNKNOWN))
        assert out["error"]["code"] == "DATASET_HANDLE_NOT_FOUND", tool.__name__
        assert out["error"]["details"]["dataset_handle"] == UNKNOWN


def test_preprocess_on_handle_leaves_global_and_other_handle_untouched(pos_mztab):
    mztab_tools.dataset_load(mztab_path=str(pos_mztab))
    a, b = _open(pos_mztab), _open(pos_mztab)
    out = json.loads(T.dataset_preprocess(normalize="pqn", dataset_handle=a))
    assert out["status"] == "success"
    assert handles.registry.get(a).dataset.pp_matrix is not None
    assert handles.registry.get(b).dataset.pp_matrix is None
    assert session_state.session.dataset.pp_matrix is None


def test_differential_and_search_on_handle(pos_mztab):
    a = _open(pos_mztab)
    T.dataset_preprocess(dataset_handle=a)
    diff = json.loads(T.dataset_differential(group_a=HC, group_b=RA, dataset_handle=a))
    assert diff["status"] == "success"
    search = json.loads(T.dataset_feature_search(feature_ids=["1"], dataset_handle=a))
    assert search["stats_source"]["result_ids"] == [diff["result_id"]]
    status = json.loads(mztab_tools.dataset_status(dataset_handle=a))
    assert status["n_features"] == 3
    assert session_state.session.dataset is None


def test_library_tools_on_handle_without_library_are_missing_state(pos_mztab):
    a = _open(pos_mztab)
    out = json.loads(library_tools.library_match_feature(precursor_mz=104.107, dataset_handle=a))
    assert out["error"]["code"] == "missing_state"


def test_library_load_into_handle_keeps_global_empty(monkeypatch, tmp_path, pos_mztab):
    msp = tmp_path / "lab" / "pos.msp"
    msp.parent.mkdir()
    msp.write_text("NAME: Choline\nPRECURSORMZ: 104.107\nIONMODE: Positive\nNum Peaks: 2\n"
                   "60.0808 100\n104.107 50\n", encoding="ascii")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(msp))
    monkeypatch.setenv("LIPIDMIX_LIBRARY_CACHE_DIR", str(tmp_path / "cache"))
    a = _open(pos_mztab)
    out = json.loads(library_tools.library_load(ion_mode="positive", dataset_handle=a))
    assert out["status"] == "success"
    assert str(msp.parent) not in json.dumps(out)
    assert handles.registry.get(a).library.store is not None
    assert session_state.session.library.store is None


def _scrub(value, tmp: str):
    """乱数の ID・一時ディレクトリ・パスを含む指紋を伏せる（それ以外は 1 文字も変えずに比べる）。"""
    import re
    text = json.dumps(value, ensure_ascii=False, sort_keys=True)
    text = re.sub(r"res_[0-9a-f]{32}", "res_X", text)
    # データセットの指紋は読んだファイルのパスを含む（置き場所で変わる）
    text = re.sub(r'"dataset_hash": "[0-9a-f]{64}"', '"dataset_hash": "X"', text)
    for form in (tmp, tmp.replace("\\", "/"), json.dumps(tmp)[1:-1]):
        text = text.replace(form, "<TMP>")
    return json.loads(text)


def golden_outputs(directory: Path) -> dict:
    """ハンドルを渡さない呼び出しの戻り値（変更前のコードでも同じ関数で採る）。"""
    path = write_mztab(directory / "AlignResult-1.mzTab")
    session_state.session = session_state.AnalysisSession()
    mztab_tools.dataset_load(mztab_path=str(path))

    def load(text, drop=()):
        data = json.loads(text)
        for key in drop:
            data.pop(key, None)
        return data

    out = {
        "status": load(mztab_tools.dataset_status(), ["server_version", "update_available"]),
        "search_window": load(T.dataset_feature_search(min_mz=100, max_mz=600)),
        "preprocess": load(T.dataset_preprocess(normalize="pqn")),
        "differential": load(T.dataset_differential(group_a=HC, group_b=RA)),
        "search_targets": load(T.dataset_feature_search(targets=[{
            "target_id": "t1", "inchikey": "OEYIOHPDSNJKLS-UHFFFAOYSA-N",
            "expected_mz": 104.107, "adduct": "[M]+", "rt_min": 1.2}])),
    }
    return _scrub(out, str(directory))


@pytest.mark.parametrize("rdkit", [True, False])
def test_without_handle_output_is_unchanged(tmp_path, monkeypatch, rdkit):
    """ハンドルを渡さない呼び出しの戻り値は、変更前のコードで採った値と同じ。

    golden は RDKit が読める機械で採った。fixture に SMILES / InChI が無いので RDKit の有無は
    `inchikey_coverage.rdkit_available` の旗にしか出ない——旗を固定し、どちらでも比べる
    （Windows のアプリケーション制御が RDKit の DLL を拒否する機械でも落ちないように）。
    """
    monkeypatch.setattr(mztab_identity, "rdkit_available", lambda: rdkit)
    golden = json.loads(GOLDEN.read_text(encoding="utf-8"))
    golden["status"]["inchikey_coverage"]["rdkit_available"] = rdkit
    assert golden_outputs(tmp_path) == golden


# ---------- スレッドとロック（mcp.call_tool 経由） ----------

def test_same_handle_calls_are_serialized(pos_mztab):
    handle_id = _open(pos_mztab)
    handle = handles.registry.get(handle_id)
    real_lock = handle.lock
    spans = []

    class SpyLock:
        def __enter__(self):
            real_lock.acquire()
            start = time.perf_counter()
            time.sleep(0.2)
            spans.append(start)
            return self

        def __exit__(self, *exc):
            spans.append(time.perf_counter())
            real_lock.release()

    handle.lock = SpyLock()

    async def main():
        await asyncio.gather(
            server.mcp.call_tool("dataset_status", {"dataset_handle": handle_id}),
            server.mcp.call_tool("dataset_status", {"dataset_handle": handle_id}))

    asyncio.run(main())
    start1, end1, start2, end2 = spans[0], spans[1], spans[2], spans[3]
    assert end1 <= start2


def test_call_without_handle_stays_on_event_loop_thread(monkeypatch):
    seen = []
    original = mztab_tools.missing_state

    def spy(*args, **kwargs):
        seen.append(threading.current_thread())
        return original(*args, **kwargs)

    monkeypatch.setattr(mztab_tools, "missing_state", spy)

    async def main():
        seen.append(threading.current_thread())
        await server.mcp.call_tool("dataset_status", {})

    asyncio.run(main())
    assert len(seen) == 2 and seen[0] is seen[1]
