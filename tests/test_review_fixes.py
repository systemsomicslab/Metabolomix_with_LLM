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
from metabolomix.core import handles, mcp_core, session_state
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


# ---------- I-2: アラインメントのファイルが無いハンドルは対話のフォルダへ落ちない ----------

@pytest.fixture
def interactive_folder(tmp_path, monkeypatch):
    """対話のデータフォルダに、別のデータセットの .arf2 / .dcl / ライブラリを置く。"""
    folder = tmp_path / "interactive"
    folder.mkdir()
    (folder / "AlignmentResult_9.arf2").write_bytes(b"placeholder")
    (folder / "AlignmentResult_9.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": 104.107, "rt": 1.2, "spectrum": [(60.08, 999.0)]},
        {"precursor_mz": 104.107, "rt": 1.2, "spectrum": [(60.08, 999.0)]},
    ]))
    monkeypatch.setattr(mcp_core, "DATA_DIR", folder)
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _p: [
        {"AlignmentID": 1, "MassCenter": 104.107, "RT": 1.2, "IonMode": "Positive"}])
    lab = tmp_path / "lab"
    lab.mkdir()
    for name, mode in (("pos.msp", "Positive"), ("neg.msp", "Negative")):
        (lab / name).write_text(f"NAME: x\nPRECURSORMZ: 104.107\nIONMODE: {mode}\n"
                                "Num Peaks: 1\n60.08 100\n", encoding="ascii")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(lab / "pos.msp"),
                                       "library.msp_negative": str(lab / "neg.msp")})
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    return folder


def _handle(tmp_path, library="auto"):
    path = write_mztab(tmp_path / "handle_data" / "AlignResult-1.mzTab")
    out = json.loads(handle_tools.load_dataset(mztab_path=str(path), library=library))
    assert out["arf2"] is None and out["dcl"] is None
    if library == "auto":
        # mzTab-M の記憶は `unknown` なので、library_load は極性の設定（[library] msp_positive）から引く
        # （mzTab-M の本番の経路）。
        library_tools.library_load(dataset_handle=out["dataset_handle"])
    return out["dataset_handle"]


def test_spot_spectrum_on_handle_without_alignment_files(tmp_path, interactive_folder):
    handle = _handle(tmp_path)
    assert json.loads(dcl_tools.dcl_spot_spectrum(spot_id=1))["status"] == "ok"   # 対話の経路
    out = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=1, dataset_handle=handle))
    assert out["status"] == "no_alignment_files"


def test_match_on_handle_without_alignment_files(tmp_path, interactive_folder):
    handle = _handle(tmp_path)
    by_spot = json.loads(library_tools.library_match_feature(spot_id=1, dataset_handle=handle))
    assert by_spot["status"] == "no_alignment_files"
    by_mz = json.loads(library_tools.library_match_feature(precursor_mz=104.107,
                                                           dataset_handle=handle))
    assert by_mz["status"] == "no_alignment_files"


def test_diagnostic_ions_on_handle_without_alignment_files(tmp_path, interactive_folder):
    handle = _handle(tmp_path)
    out = json.loads(dcl_tools.dcl_diagnostic_ions([1], [{"mz": 60.08}], dataset_handle=handle))
    assert "\tno_alignment_files" in out["ions"]


def test_library_load_on_handle_uses_handle_polarity(tmp_path, interactive_folder):
    handle = _handle(tmp_path, library="none")
    out = json.loads(library_tools.library_load(dataset_handle=handle))
    assert out["status"] == "success" and out["file"] == "pos.msp"
    wrong = json.loads(library_tools.library_load(ion_mode="negative", dataset_handle=handle))
    assert wrong["code"] == "LIBRARY_POLARITY_MISMATCH"


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
