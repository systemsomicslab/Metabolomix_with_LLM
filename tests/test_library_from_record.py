"""下流の照合は、データセットが覚えたライブラリを既定にする（B-1 Task 10）。

ライブラリのファイル名は作り物（研究室の実ファイルに依存しない・書かない）。
"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import handles, mcp_core, session_state
from metabolomix.library import tools as library_tools
from metabolomix.tools import handle_tools
from metabolomix.tools import reanalysis_tools as R
from tests.config_fixtures import use_config
from tests.dbs_fixture import record, write_dbs
from tests.dcl_fixture import build_dcl_bytes
from tests.msdial_folder_fixtures import write_msdial_folder

_SAMPLES = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(3)]
_STEM = "AlignmentResult_2026_02_01_00_00_00"


@pytest.fixture(autouse=True)
def fresh_state(tmp_path, monkeypatch):
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path / "datadir")
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _folder(tmp_path, name="gui"):
    folder = tmp_path / name
    write_msdial_folder(folder, spots=[{"name": "GABA", "mz": 500.0, "rt": 1.23}],
                        samples=_SAMPLES)
    return folder


def _dbs(path, mz=500.0):
    return write_dbs(path, [("lib", [record(mz=mz, ion_mode=0)], ["k"])])


def _open(folder, **kwargs) -> dict:
    return json.loads(handle_tools.dataset_handle_open(directory=str(folder), **kwargs))


def test_handle_opens_remembered_library(tmp_path):
    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    out = _open(folder)
    assert out["status"] == "success", out
    assert out["library"]["file"] == "Proj_Loaded.msp2.dbs" and out["library"]["record_count"] == 1
    assert "library_note" not in out
    assert str(folder) not in json.dumps(out, ensure_ascii=False)
    store = handles.registry.get(out["dataset_handle"]).library.store
    with pytest.raises(Exception):  # 読み取り専用
        store._conn.execute("INSERT INTO meta(key, value) VALUES ('x', '1')")


def test_library_sha_mismatch_not_opened(tmp_path):
    folder = _folder(tmp_path)
    dbs = _dbs(folder / "Proj_Loaded.msp2.dbs")
    first = _open(folder, library="none")  # 記憶（sha256）はこの時点で測られる
    ds = handles.registry.get(first["dataset_handle"]).dataset
    assert ds.library_record["status"] == "known"
    _dbs(dbs, mz=501.0)  # 開いた後でライブラリを書き換える
    store, note = handle_tools._open_library_for(ds)
    assert store is None and "不一致" in note
    assert str(folder) not in note


def test_no_record_does_not_use_polarity_setting(tmp_path, monkeypatch):
    from tests.reanalysis_fixtures import write_mztab

    lab = tmp_path / "lab"
    lab.mkdir()
    msp = lab / "pos.msp"
    msp.write_text("NAME: x\nPRECURSORMZ: 104.107\nIONMODE: Positive\nNum Peaks: 1\n60.08 100\n",
                   encoding="ascii")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(msp)})
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab")
    out = json.loads(handle_tools.dataset_handle_open(mztab_path=str(mztab)))
    assert out["polarity"] == "positive" and out["library"] is None
    assert out["library_record"]["status"] == "unknown"
    assert "極性の設定" in out["library_note"]
    assert handles.registry.get(out["dataset_handle"]).library.store is None


def test_library_load_without_args_uses_record(tmp_path):
    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    h = _open(folder, library="none")["dataset_handle"]
    out = json.loads(library_tools.library_load(dataset_handle=h))
    assert out["status"] == "success" and out["file"] == "Proj_Loaded.msp2.dbs"
    assert out["library_source"] == "dataset_record" and "differs_from_record" not in out
    assert handles.registry.get(h).library.store is not None
    assert str(folder) not in json.dumps(out, ensure_ascii=False)


def test_library_load_without_args_uses_latest_dataset_record(tmp_path):
    from metabolomix.tools import mztab_tools

    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    mztab_tools.dataset_load(directory=str(folder))
    assert session_state.session.dataset is not None
    out = json.loads(library_tools.library_load())
    assert out["status"] == "success" and out["library_source"] == "dataset_record"
    assert out["file"] == "Proj_Loaded.msp2.dbs"
    assert session_state.session.library.store is not None


def test_explicit_file_path_overrides_and_is_reported(tmp_path):
    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    (tmp_path / "other").mkdir()
    other = _dbs(tmp_path / "other" / "Other_Loaded.msp2.dbs", mz=600.0)
    h = _open(folder, library="none")["dataset_handle"]
    out = json.loads(library_tools.library_load(file_path=str(other), dataset_handle=h))
    assert out["status"] == "success" and out["file"] == "Other_Loaded.msp2.dbs"
    assert out["differs_from_record"] == {"remembered": "Proj_Loaded.msp2.dbs",
                                          "opened": "Other_Loaded.msp2.dbs"}
    assert str(other.parent) not in json.dumps(out, ensure_ascii=False)


def test_explicit_same_file_is_not_reported_as_different(tmp_path):
    folder = _folder(tmp_path)
    dbs = _dbs(folder / "Proj_Loaded.msp2.dbs")
    h = _open(folder, library="none")["dataset_handle"]
    out = json.loads(library_tools.library_load(file_path=str(dbs), dataset_handle=h))
    assert out["status"] == "success" and "differs_from_record" not in out


def test_library_load_mismatch_does_not_fall_back(tmp_path):
    folder = _folder(tmp_path)
    dbs = _dbs(folder / "Proj_Loaded.msp2.dbs")
    h = _open(folder, library="none")["dataset_handle"]
    _dbs(dbs, mz=501.0)
    out = json.loads(library_tools.library_load(dataset_handle=h))
    assert out["status"] == "error" and "不一致" in out["message"] and "file_path" in out["message"]
    assert handles.registry.get(h).library.store is None


def test_library_load_lbm_only_does_not_fall_back(tmp_path):
    folder = _folder(tmp_path)
    write_dbs(folder / "Proj_Loaded.lbm2.dbs", [("lib", [record()], ["k"])])
    h = _open(folder)["dataset_handle"]
    out = json.loads(library_tools.library_load(dataset_handle=h))
    assert out["status"] == "error" and "LBM" in out["message"] and "file_path" in out["message"]


def test_lbm_only_dataset_returns_not_evaluable(tmp_path):
    folder = _folder(tmp_path)
    write_dbs(folder / "Proj_Loaded.lbm2.dbs", [("lib", [record()], ["k"])])
    opened = _open(folder)
    assert opened["library"] is None and opened["library_record"]["status"] == "lbm_only"
    h = opened["dataset_handle"]
    match = json.loads(library_tools.library_match_feature(precursor_mz=500.0, dataset_handle=h))
    assert match["status"] == "not_evaluable" and match["reason"]
    assert "error" not in match
    qualify = json.loads(R.library_msms_qualify({"target": {"inchikey": "X"}}, precursor_mz=500.0,
                                                dataset_handle=h))
    assert qualify["status"] == "not_evaluable" and qualify["verdict"] == "not_evaluable"
    assert str(folder) not in json.dumps([match, qualify], ensure_ascii=False)


def test_no_library_and_no_lbm_record_is_still_missing_state(tmp_path):
    folder = _folder(tmp_path)
    h = _open(folder)["dataset_handle"]
    out = json.loads(library_tools.library_match_feature(precursor_mz=500.0, dataset_handle=h))
    assert out["error"]["code"] == "missing_state"


def _matching_folder(tmp_path, *, with_lbm: bool):
    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    if with_lbm:
        write_dbs(folder / "Proj_Loaded.lbm2.dbs", [("lib", [record()], ["k"])])
    entries = [{"precursor_mz": 500.0, "rt": 1.23,
                "spectrum": [(87.04, 999.0), (69.03, 500.0), (500.0, 100.0)]}]
    dcl = folder / f"{_STEM}.dcl"
    dcl.write_bytes(build_dcl_bytes(entries))
    return folder, dcl


@pytest.mark.parametrize("with_lbm, annotator, expected", [
    (True, "C:/lab/Proj.lbm2", True),     # LBM も使い、代表の注釈器が LBM
    (True, "LabMsp", False),              # LBM も使ったが、代表は LBM ではない
    (False, "C:/lab/Proj.lbm2", False),   # LBM は使っていない
])
def test_identified_by_lbm_only_when_lbm_also_and_annotator_is_lbm(tmp_path, with_lbm,
                                                                    annotator, expected):
    folder, dcl = _matching_folder(tmp_path, with_lbm=with_lbm)
    opened = _open(folder)
    handle = handles.registry.get(opened["dataset_handle"])
    assert handle.library.store is not None
    fid = handle.dataset.feature_ids[0]
    handle.dataset.feature_metadata[fid]["annotator_id"] = annotator
    out = json.loads(library_tools.library_match_feature(
        spot_id=0, dcl_file=str(dcl),
        dataset_handle=opened["dataset_handle"]))
    assert out["status"] == "success", out
    if expected:
        assert out["identified_by_lbm"] is True
    else:
        assert "identified_by_lbm" not in out  # False も書かない


# ---------- fix round 1 ----------

def _msp(path, name="x"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"NAME: {name}\nPRECURSORMZ: 500.0\nIONMODE: Positive\nNum Peaks: 1\n60.08 100\n",
                    encoding="ascii")
    return path


def test_pipeline_pinned_sha_differs_from_current_file_is_not_opened(tmp_path):
    """Review Focus 4: run が固定した sha256 と今のファイルが違えば、開かず「不一致」。"""
    from tests.test_dataset_opening import _job, _pipeline

    lib = _msp(tmp_path / "libs" / "Run_Lab.msp")
    run = _pipeline(tmp_path / "p", _job(tmp_path / "p"))
    record_path = run / "pipeline-run.json"
    data = json.loads(record_path.read_text(encoding="utf-8"))
    data["inputs"] = {"library_mode": "msp_only", "msp": {
        "path": str(lib), "file": lib.name, "source": "config_file", "sha256": "0" * 64}}
    record_path.write_text(json.dumps(data), encoding="utf-8")
    out = json.loads(handle_tools.dataset_handle_open(pipeline_path=str(run), library="auto"))
    assert out["status"] == "success", out
    assert out["library"] is None and "不一致" in out["library_note"]
    assert out["library_record"]["status"] == "mismatch"
    assert handles.registry.get(out["dataset_handle"]).library.store is None
    assert str(lib.parent) not in json.dumps(out, ensure_ascii=False)


def test_stale_dataset_library_load_refuses_after_the_file_changes(tmp_path):
    """データセットを開いた後で記憶のファイルが変わったら、引数なしの library_load は開かない。"""
    from metabolomix.tools import mztab_tools

    folder = _folder(tmp_path)
    dbs = _dbs(folder / "Proj_Loaded.msp2.dbs")
    mztab_tools.dataset_load(directory=str(folder))
    _dbs(dbs, mz=501.0)
    out = json.loads(library_tools.library_load())
    assert out["status"] == "error" and "不一致" in out["message"]
    assert session_state.session.library.store is None


def test_ion_mode_skips_the_record_and_reports_the_difference(tmp_path, monkeypatch):
    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    setting = _msp(tmp_path / "lab" / "Setting_Pos.msp")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(setting)})
    h = _open(folder, library="none")["dataset_handle"]
    out = json.loads(library_tools.library_load(ion_mode="positive", dataset_handle=h))
    assert out["status"] == "success" and out["file"] == "Setting_Pos.msp"
    assert "library_source" not in out
    assert out["differs_from_record"] == {"remembered": "Proj_Loaded.msp2.dbs",
                                          "opened": "Setting_Pos.msp"}
    assert str(setting.parent) not in json.dumps(out, ensure_ascii=False)


def test_ion_mode_opposite_to_the_handle_polarity_still_mismatches(tmp_path, monkeypatch):
    folder = _folder(tmp_path)
    _dbs(folder / "Proj_Loaded.msp2.dbs")
    setting = _msp(tmp_path / "lab" / "Setting_Neg.msp")
    use_config(monkeypatch, tmp_path, {"library.msp_negative": str(setting)})
    h = _open(folder, library="none")["dataset_handle"]
    handles.registry.get(h).polarity = "positive"
    out = json.loads(library_tools.library_load(ion_mode="negative", dataset_handle=h))
    assert out["code"] == "LIBRARY_POLARITY_MISMATCH"


def test_identical_copy_is_not_different_but_a_rewritten_file_is(tmp_path):
    import shutil

    folder = _folder(tmp_path)
    dbs = _dbs(folder / "Proj_Loaded.msp2.dbs")
    copy = tmp_path / "elsewhere" / "Copy_Of_Loaded.msp2.dbs"
    copy.parent.mkdir()
    shutil.copy(dbs, copy)
    h = _open(folder, library="none")["dataset_handle"]
    same = json.loads(library_tools.library_load(file_path=str(copy), dataset_handle=h))
    assert same["status"] == "success" and "differs_from_record" not in same
    _dbs(dbs, mz=501.0)  # 同じ場所で中身が変わった
    changed = json.loads(library_tools.library_load(file_path=str(dbs), dataset_handle=h))
    assert changed["status"] == "success"
    assert changed["differs_from_record"] == {"remembered": "Proj_Loaded.msp2.dbs",
                                              "opened": "Proj_Loaded.msp2.dbs"}


def test_lbm_only_with_explicit_file_path_says_it_was_not_identified_with_it(tmp_path):
    folder = _folder(tmp_path)
    write_dbs(folder / "Proj_Loaded.lbm2.dbs", [("lib", [record()], ["k"])])
    (tmp_path / "other").mkdir()
    other = _dbs(tmp_path / "other" / "Other.msp2.dbs")
    h = _open(folder)["dataset_handle"]
    out = json.loads(library_tools.library_load(file_path=str(other), dataset_handle=h))
    assert out["status"] == "success" and "LBM" in out["note"] and "Other.msp2.dbs" in out["note"]
    assert "differs_from_record" not in out


def test_session_dataset_without_a_record_falls_back_to_the_settings(tmp_path, monkeypatch):
    """R21: mzTab-M の記憶は unknown。引数なしの library_load は極性の設定へ落ちる。"""
    from metabolomix.tools import mztab_tools
    from tests.reanalysis_fixtures import write_mztab

    setting = _msp(tmp_path / "lab" / "Setting_Pos.msp")
    use_config(monkeypatch, tmp_path, {"library.msp_positive": str(setting)})
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab")
    mztab_tools.dataset_load(mztab_path=str(mztab))
    assert session_state.session.dataset.library_record["status"] == "unknown"
    out = json.loads(library_tools.library_load())
    assert out["status"] == "success" and out["file"] == "Setting_Pos.msp"
    assert "library_source" not in out
