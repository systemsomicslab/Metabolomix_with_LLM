"""全体レビューの指摘（I-1〜I-5）を固定する。"""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from metabolomix.analysis import analysis_spec as A
from metabolomix.analysis.recipe_matrix import knn_impute
from metabolomix.arf2 import reader as arf2_reader
from metabolomix.console import method_file as method_file_mod
from metabolomix.core import handles, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.dcl import tools as dcl_tools
from metabolomix.library import tools as library_tools
from metabolomix.msdial import alignment_spot
from metabolomix.tools import handle_tools
from tests.dcl_fixture import build_dcl_bytes
from tests.reanalysis_fixtures import write_mztab
from tests.config_fixtures import use_config

NAN = math.nan


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


# ---------- I-1: 外部の行列の重複した列 ----------

def test_external_matrix_duplicate_sample_columns_rejected(tmp_path):
    matrix = tmp_path / "m.tsv"
    matrix.write_text("feature_id\tA1\tA2\tB1\tB2\tB2\nF1\t1\t2\t3\t4\t5\n", encoding="utf-8")
    groups = tmp_path / "g.tsv"
    groups.write_text("sample\tgroup\nA1\tA\nA2\tA\nB1\tB\nB2\tB\n", encoding="utf-8")
    spec = {"input": {"source": "tsv", "matrix_path": str(matrix)},
            "groups": {"source": "tsv", "path": str(groups)},
            "comparison": {"reference": "A", "test": "B"},
            "analyses": [{"analysis_id": "m", "test": {"kind": "welch"}}]}
    with pytest.raises(DomainError) as info:
        A.run_spec(A.normalize_spec(spec))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"
    assert info.value.details["duplicated"] == ["B2"]


def test_external_matrix_empty_id_rejected(tmp_path):
    matrix = tmp_path / "m.tsv"
    matrix.write_text("feature_id\tA1\tA2\n\t1\t2\nF2\t3\t4\n", encoding="utf-8")
    groups = tmp_path / "g.tsv"
    groups.write_text("sample\tgroup\nA1\tA\nA2\tA\n", encoding="utf-8")
    spec = {"input": {"source": "tsv", "matrix_path": str(matrix)},
            "groups": {"source": "tsv", "path": str(groups)},
            "comparison": {"reference": "A", "test": "B"},
            "analyses": [{"analysis_id": "m", "test": {"kind": "welch"}}]}
    with pytest.raises(DomainError) as info:
        A.run_spec(A.normalize_spec(spec))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"


# ---------- I-2: アラインメントのファイルが無いハンドルは、別のデータセットのファイルへ落ちない ----------

@pytest.fixture
def interactive_folder(tmp_path, monkeypatch):
    """別のデータセット（アラインメントの .arf2 / .dcl とライブラリを持つ解析フォルダ）を 1 つ開いておく。

    ハンドルを渡さない呼び出し・ほかのハンドルの呼び出しが、このデータセットのファイルへ落ちないことを
    確かめるための「落ちうる先」。ハンドルの ID を `handle` に持つ。"""
    from tests.curation_fixtures import write_alignment_set

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    paths = write_alignment_set(tmp_path / "interactive")
    out = json.loads(handle_tools.load_dataset(directory=str(paths["arf2"].parent), library="none"))
    assert out["status"] == "success", out
    return {**paths, "handle": out["dataset_handle"]}


def _handle(tmp_path, interactive_folder):
    """`.arf2` / `.dcl` を持たない（mzTab-M だけの）ハンドル。ライブラリは `file_path` で開く。
    直近になる（`load_dataset` はいつも直近にする）。"""
    path = write_mztab(tmp_path / "handle_data" / "AlignResult-1.mzTab")
    out = json.loads(handle_tools.load_dataset(mztab_path=str(path), library="none"))
    assert out["arf2"] is None and out["dcl"] is None
    loaded = json.loads(library_tools.library_load(file_path=str(interactive_folder["msp"]),
                                                   dataset_handle=out["dataset_handle"]))
    assert loaded["status"] == "success", loaded
    return out["dataset_handle"]


def test_spot_spectrum_on_handle_without_alignment_files(tmp_path, interactive_folder):
    other = interactive_folder["handle"]
    assert json.loads(dcl_tools.dcl_spot_spectrum(spot_id=0, dataset_handle=other))["status"] == "ok"
    handle = _handle(tmp_path, interactive_folder)
    out = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=0, dataset_handle=handle))
    assert out["status"] == "no_alignment_files"
    # ハンドルを省いた呼び出し（直近はファイルを持たないハンドル）も、開いている別のハンドルへは落ちない。
    assert json.loads(dcl_tools.dcl_spot_spectrum(spot_id=0))["status"] == "no_alignment_files"


def test_match_on_handle_without_alignment_files(tmp_path, interactive_folder):
    handle = _handle(tmp_path, interactive_folder)
    by_spot = json.loads(library_tools.library_match_feature(spot_id=0, dataset_handle=handle))
    assert by_spot["status"] == "no_alignment_files"
    by_mz = json.loads(library_tools.library_match_feature(precursor_mz=760.5851,
                                                           dataset_handle=handle))
    assert by_mz["status"] == "no_alignment_files"


def test_diagnostic_ions_on_handle_without_alignment_files(tmp_path, interactive_folder):
    handle = _handle(tmp_path, interactive_folder)
    out = json.loads(dcl_tools.dcl_diagnostic_ions([0], [{"mz": 184.07}], dataset_handle=handle))
    assert "	no_alignment_files" in out["ions"]


# ---------- I-3: .arf2 のカタログのキャッシュ ----------

def test_catalog_cache_keeps_several_files(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(arf2_reader, "deserialize", lambda stream: calls.append(1) or [])
    a, b = tmp_path / "a.arf2", tmp_path / "b.arf2"
    a.write_bytes(b"a")
    b.write_bytes(b"b")
    for path in (a, b, a, b, a):
        arf2_reader.load_catalog(path)
    assert len(calls) == 2


# ---------- I-4: KNN の欠測の多い行（R impute.knn の knnimp と同じ扱い） ----------

def test_knn_heavy_rows_are_not_neighbours_and_use_heavy_block_means():
    data = np.array([
        [1, 2, 3, NAN],          # 軽い行（欠測 1/4）
        [1, NAN, NAN, 100],      # 重い行（2/4 > rowmax 0.4）——近傍にしない
        [NAN, 7, 5, NAN],        # 重い行
        [1, 2, 3, 4],            # 軽い行（完全）
    ], dtype=float)
    out = knn_impute(data, k=1, rowmax=0.4, colmax=1.0)
    assert out[0, 3] == pytest.approx(4.0)                 # 重い行の 100 を借りない
    assert out[1, 1] == pytest.approx(7.0) and out[1, 2] == pytest.approx(5.0)
    assert out[2, 0] == pytest.approx(1.0) and out[2, 3] == pytest.approx(100.0)


def test_knn_fallback_column_mean_is_over_light_rows():
    data = np.array([
        [1, 2, NAN],
        [1, 2, NAN],             # 近傍も同じ列が欠測
        [50, 60, 70],            # 遠い軽い行
        [NAN, NAN, 1000],        # 重い行（2/3 > 0.5）——列の平均に入れない
    ], dtype=float)
    out = knn_impute(data, k=1, rowmax=0.5, colmax=1.0)
    assert out[0, 2] == pytest.approx(70.0)                # 軽い行の列 2 の観測値の平均


# ---------- I-5: コメント行の ASCII 以外 ----------

def test_comment_lines_with_non_ascii_are_allowed(tmp_path):
    src = tmp_path / "method.txt"
    src.write_bytes("Ion mode: Positive\n# 日本語のメモ\nMS1 tolerance: 0.01\n".encode("utf-8"))
    method_file_mod.check_ascii(src)                       # コメントは Console が読まない
    out = method_file_mod.write_effective_method_file(src, tmp_path / "out.txt",
                                                      {"Ion mode": "Negative"})
    lines = out.read_text(encoding="ascii").splitlines()
    assert lines[0] == "Ion mode: Negative" and lines[1].startswith("# ")
