"""入口 `open_dataset`（解析フォルダ・mzTab-M・job・pipeline を 1 つの入口で開く）と、
既存の入口ツール 2 つ（`load_dataset` / `load_dataset`）への配線（B-1 Task 8）。

fixture はすべてこのテストが tmp に作る（研究室のライブラリ・研究データに依存しない）。
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from metabolomix.core import handles, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.opening import open_dataset
from metabolomix.handoff.schema import SCHEMA_VERSION
from metabolomix.mztab.loading import load_dataset_state
from tests.curation_fixtures import match_result
from tests.dbs_fixture import record, write_dbs
from tests.msdial_folder_fixtures import write_matching_mztab, write_msdial_folder

_STEM = "AlignmentResult_2026_02_01_00_00_00"


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _spec() -> dict:
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
               for i in range(4)]
    spots = [
        {"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0],
         "gap_filled": [False, True, False, False], "msms": [True, False, False, False],
         "matches": [match_result({0: "PC 34:1"})]},
        {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]},
    ]
    return {"spots": spots, "samples": samples}


def _folder(path: Path, **kwargs) -> Path:
    write_msdial_folder(path, **_spec(), **kwargs)
    return path


def _job(tmp_path: Path, *, status: str = "completed", with_folder: bool = True,
         method_file: str = "m.txt", spec: dict | None = None,
         with_mddata: bool = True, mztab_name: str = "Height_x.mzTab") -> Path:
    """Console の run の形（run_dir/msdial/ に mzTab-M と、あれば `.arf` 一式）の job。"""
    spec = spec or _spec()
    run_dir = tmp_path / "runs" / "job1"
    out = run_dir / "msdial"
    if with_folder:
        write_msdial_folder(out, **spec, with_mddata=with_mddata)
    write_matching_mztab(out / mztab_name, spec)
    job_path = run_dir / "analysis-job.json"
    job_path.write_text(json.dumps({
        "schema": SCHEMA_VERSION, "job_id": "job1", "status": status,
        "created_at": "t", "updated_at": "t",
        "source": {"dataset_root": str(tmp_path), "input_count": 4},
        "software": {"name": "MS-DIAL", "version": "5", "execution_mode": "console",
                     "method_file": method_file},
        "project": {"omics": "lipidomics", "polarity": "positive", "measure": "peak_height"},
        "run_dir": str(run_dir),
        "primary_mztab_files": [{"path": f"msdial/{mztab_name}", "polarity": "positive",
                                 "measure": "peak_height", "sha256": "", "validation": {},
                                 "root": "run_dir"}],
        "artifacts": [], "sample_manifest": None, "warnings": [], "error": None,
    }), encoding="utf-8")
    return job_path


def _manifest_row(name: str, group: str) -> dict:
    """run が解決した試料対応表の 1 行（`.mddata` の raw のパスで試料に結ぶ）。"""
    values = {"sample_id": name, "source_file": f"C:/data/{name}.wiff", "role": "sample",
              "group": group, "batch": "1", "injection_order": None, "qc_pool": None,
              "include": True}
    declared = {"source": "manifest", "confidence": "declared"}
    provenance = {field: {"value": value, **declared} for field, value in values.items()}
    return {**values, "provenance": provenance, "conflicts": []}


def _pipeline(tmp_path: Path, job_path: Path, *, rows: list[dict] | None = None) -> Path:
    root = tmp_path / "pipeline_x"
    root.mkdir()
    inputs = {"library_mode": "auto",
              "msp": {"path": None, "file": None, "sha256": None, "source": "not_used"}}
    if rows is not None:
        inputs["manifest"] = rows
    (root / "pipeline-run.json").write_text(json.dumps({
        "schema": "pipeline-run.v1", "inputs": inputs,
        "upstream": {"console_job_path": str(job_path)}, "results": []}), encoding="utf-8")
    return root


def _code(call) -> DomainError:
    with pytest.raises(DomainError) as exc:
        call()
    return exc.value


# ---------- 入口の選び方 ----------

def test_directory_opens_folder_adapter(tmp_path):
    folder = _folder(tmp_path / "gui")
    opened = open_dataset(directory=str(folder))
    assert opened.adapter == "msdial_folder" and opened.kind == "directory"
    assert opened.dataset.source_format == "msdial_folder"
    assert opened.selection.arf.name == f"{_STEM}_PeakProperties.arf"
    assert opened.path == str(folder.resolve())


def test_exactly_one_entry_required(tmp_path):
    folder = _folder(tmp_path / "gui")
    assert _code(lambda: open_dataset()).code == "DATASET_BAD_REQUEST"
    two = _code(lambda: open_dataset(directory=str(folder), mztab_path=str(folder / "x.mzTab")))
    assert two.code == "DATASET_BAD_REQUEST"


def test_source_conflicts_rejected(tmp_path):
    folder = _folder(tmp_path / "gui")
    mztab = write_matching_mztab(tmp_path / "direct" / "x.mzTab", _spec())
    assert _code(lambda: open_dataset(directory=str(folder), source="mztab")).code \
        == "DATASET_BAD_REQUEST"
    assert _code(lambda: open_dataset(mztab_path=str(mztab), source="folder")).code \
        == "DATASET_BAD_REQUEST"
    assert _code(lambda: open_dataset(directory=str(folder), source="arf")).code \
        == "DATASET_BAD_REQUEST"
    # measure は解析フォルダのアダプタだけが使う（mzTab-M の定量種別はファイルで決まる）。
    assert _code(lambda: open_dataset(mztab_path=str(mztab),
                                      measure="peak_area_above_zero")).code \
        == "DATASET_BAD_REQUEST"


def test_job_defaults_to_folder_and_mztab_on_request(tmp_path):
    job = _job(tmp_path)
    folder = open_dataset(job_path=str(job))
    mztab = open_dataset(job_path=str(job), source="mztab")
    assert folder.adapter == "msdial_folder"
    assert mztab.adapter == "mztab"
    assert open_dataset(job_path=str(job), source="folder").adapter == "msdial_folder"
    # 門の結果（job_path・出所の裏取り・探索専用）は解析フォルダの読み込みにも写す。
    for opened in (folder, mztab):
        assert opened.dataset.job_path == str(job)
        assert opened.dataset.source_verification == "legacy_unverified"
        assert opened.dataset.exploratory_only is False
    assert folder.selection.arf.parent == job.parent / "msdial"


def test_incomplete_job_gate_shared(tmp_path):
    job = _job(tmp_path, status="partial")
    for source in ("auto", "mztab"):
        err = _code(lambda: open_dataset(job_path=str(job), source=source))
        assert err.code == "INCOMPLETE_ANALYSIS_JOB", source
    for source in ("auto", "mztab"):
        opened = open_dataset(job_path=str(job), source=source, allow_incomplete=True)
        assert opened.dataset.exploratory_only is True, source


def test_job_without_arf_does_not_fall_back_to_mztab(tmp_path):
    job = _job(tmp_path, with_folder=False)
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert err.code == "ALIGNMENT_NOT_FOUND"
    assert 'source="mztab"' in err.message
    assert open_dataset(job_path=str(job), source="mztab").adapter == "mztab"


def test_mztab_source_matches_current_loader(tmp_path):
    from metabolomix.pipeline import session_handoff

    job = _job(tmp_path)
    run = _pipeline(tmp_path, job)
    mztab = job.parent / "msdial" / "Height_x.mzTab"
    pairs = [
        (open_dataset(job_path=str(job), source="mztab").dataset,
         load_dataset_state(job_path=job)),
        (open_dataset(mztab_path=str(mztab)).dataset, load_dataset_state(mztab_path=mztab)),
        (open_dataset(pipeline_path=str(run), source="mztab").dataset,
         session_handoff.load_run_dataset(str(run))),
    ]
    for opened, current in pairs:
        np.testing.assert_array_equal(opened.feature_matrix, current.feature_matrix)
        assert opened.feature_ids == current.feature_ids
        assert opened.sample_names == current.sample_names
        assert opened.source_format == current.source_format == "mztab"
        assert (opened.detected_mask is None) == (current.detected_mask is None)
        if current.detected_mask is not None:
            np.testing.assert_array_equal(opened.detected_mask, current.detected_mask)
        assert opened.validation_result.get("warnings") == current.validation_result.get("warnings")
        assert opened.source_verification == current.source_verification
        assert opened.exploratory_only == current.exploratory_only


def test_base_cache_returns_independent_copies(tmp_path):
    folder = _folder(tmp_path / "gui")
    first = open_dataset(directory=str(folder))
    second = open_dataset(directory=str(folder))
    assert first.reused_base is False and second.reused_base is True
    assert first.dataset is not second.dataset
    assert first.dataset.feature_matrix is not second.dataset.feature_matrix
    assert first.dataset.dataset_id != second.dataset.dataset_id
    first.dataset.feature_matrix[0, 0] = -1.0
    first.dataset.results["x"] = {}
    third = open_dataset(directory=str(folder))
    assert third.reused_base is True
    assert third.dataset.feature_matrix[0, 0] == 100.0 and third.dataset.results == {}


def test_base_cache_key_carries_allow_incomplete(tmp_path):
    """allow_incomplete=True で読んだ基を、指定しない呼び出しへ渡さない（門を迂回させない）。"""
    job = _job(tmp_path, status="partial")
    for source in ("auto", "mztab"):
        assert open_dataset(job_path=str(job), source=source,
                            allow_incomplete=True).dataset.exploratory_only is True
        err = _code(lambda: open_dataset(job_path=str(job), source=source))
        assert err.code == "INCOMPLETE_ANALYSIS_JOB", source


@pytest.mark.parametrize("sidecar", ["curation/flags.jsonl", f"{_STEM}_tags.xml"])
def test_base_cache_reloads_when_a_folder_sidecar_changes(tmp_path, sidecar):
    """解析フォルダのアダプタはタグとキュレーションのフラグも読むので、それが変われば読み直す。

    既にあるファイルの中身だけを変える（フォルダの更新時刻は変わらない）——入力のフォルダの
    指紋ではなく、サイドカー自体を指紋に数えていることを見る。
    """
    folder = _folder(tmp_path / "gui")
    path = folder / sidecar
    path.parent.mkdir(exist_ok=True)
    path.write_text("", encoding="utf-8")
    open_dataset(directory=str(folder))
    assert open_dataset(directory=str(folder)).reused_base is True
    with open(path, "a", encoding="utf-8") as handle:
        handle.write("\n")
    assert open_dataset(directory=str(folder)).reused_base is False


# ---------- pipeline ----------

def test_pipeline_folder_reads_run_manifest_but_not_matrices(tmp_path):
    job = _job(tmp_path)
    rows = [_manifest_row(f"s{i}", "ctrl" if i < 2 else "case") for i in range(4)]
    run = _pipeline(tmp_path, job, rows=rows)
    opened = open_dataset(pipeline_path=str(run))
    ds = opened.dataset
    assert opened.adapter == "msdial_folder" and opened.kind == "pipeline_path"
    assert [r["group"] for r in ds.sample_metadata_rows] == ["ctrl", "ctrl", "case", "case"]
    assert ds.analysis_matrices == {} and ds.handoff_report is None
    assert any('source="mztab"' in note for note in opened.notes)


def test_pipeline_folder_without_mddata_binds_manifest_by_mztab_raw_paths(tmp_path):
    # Console は save_project=False なので run の input/ に .mddata は無い。試料の raw のパスは
    # 正準 mzTab-M の ms_run[N]-location から与える。
    job = _job(tmp_path, with_folder=False)
    rows = [_manifest_row(f"s{i}", "ctrl" if i < 2 else "case") for i in range(4)]
    run = _pipeline(tmp_path, job, rows=rows)
    write_msdial_folder(run / "input", **_spec(), stem="AlignResult-1", with_mddata=False)
    opened = open_dataset(pipeline_path=str(run), source="folder")
    ds = opened.dataset
    assert opened.adapter == "msdial_folder" and ds.msdial_samples is None
    assert [p.replace("\\", "/") for p in ds.sample_raw_paths] == [
        f"C:/data/s{i}.wiff" for i in range(4)]
    assert [r["group"] for r in ds.sample_metadata_rows] == ["ctrl", "ctrl", "case", "case"]
    assert [r["sample_id"] for r in ds.sample_metadata_rows] == ds.sample_names
    # directory= で開くときは変わらない（job の mzTab-M が無いので raw のパスは None のまま）。
    direct = open_dataset(directory=str(run / "input")).dataset
    assert direct.sample_raw_paths is None


def _shifted_spec() -> dict:
    """50 スポット × 4 試料。"""
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
               for i in range(4)]
    spots = [{"name": "Unknown", "mz": 300.0 + 10 * k, "heights": [100.0, 110.0, 120.0, 130.0]}
             for k in range(50)]
    return {"spots": spots, "samples": samples}


def _shift_sample0_mz(arf: Path, spec: dict, spot_indices=(10, 20), delta=0.023) -> None:
    """同じアラインメントのまま、試料 0 自身のピークの m/z だけを `delta` ずらす（スポットの
    m/z = 検出セルの高さ最大のもの = 試料 3 は変わらない）。"""
    from tests.curation_fixtures import arf_row, write_arf

    groups = []
    for k, spot in enumerate(spec["spots"]):
        rows = []
        for sample, height in zip(spec["samples"], spot["heights"]):
            mz = spot["mz"] + (delta if k in spot_indices and sample["file_id"] == 0 else 0.0)
            row = arf_row(file_id=sample["file_id"], mz=mz, rt=12.0, height=height,
                          gap_filled=False)
            row[1] = sample["name"]
            rows.append(row)
        groups.append(rows)
    write_arf(arf, groups)


def test_arf_spot_mz_is_the_highest_detected_cell_not_sample_zero(tmp_path):
    spec = _shifted_spec()
    paths = write_msdial_folder(tmp_path / "gui", **spec)
    _shift_sample0_mz(paths["arf"], spec)
    # 試料 0 の行の m/z だけ 0.023 ずれたスポットが 4%。`.arf` から作ったスポットの代表 m/z
    # （検出されたセルのうち高さ最大＝試料 3）はずれていない。試料 0 の m/z（400.023）ではない。
    from metabolomix.mztab import evidence

    scan = evidence.scan_arf(paths["arf"], measures=())
    assert scan.spot_mz[10] == pytest.approx(400.0)
    assert open_dataset(directory=str(tmp_path / "gui")).dataset.feature_metadata["10"]["mz"] == 400.0


def test_pipeline_opens_when_sample_zero_mz_is_shifted(tmp_path):
    spec = _shifted_spec()
    job = _job(tmp_path, with_folder=False)
    write_matching_mztab(job.parent / "msdial" / "Height_x.mzTab", spec)
    run = _pipeline(tmp_path, job)
    paths = write_msdial_folder(run / "input", **spec, stem="AlignResult-1")
    _shift_sample0_mz(paths["arf"], spec)
    opened = open_dataset(pipeline_path=str(run), source="folder")
    assert opened.dataset.source_format == "msdial_folder"
    assert len(opened.dataset.feature_ids) == 50


def test_pipeline_input_fallback_already_a_candidate_is_listed_once(tmp_path):
    # mzTab-M の隣の `.arf` が候補にあり、run の input/ の `.arf2` も同じ `.arf` を指すとき、
    # 二重に検証・列挙しない（解決後のパスで比べる）。
    job = _job(tmp_path, with_folder=False)
    run = _pipeline(tmp_path, job)
    paths = write_msdial_folder(run / "input", **_other_spec(), stem="AlignResult-1")
    artifacts = json.loads(job.read_text(encoding="utf-8"))
    artifacts["artifacts"] = [{"role": "peak_matrix_source", "path": str(paths["arf"]),
                               "format": "arf", "root": "run_dir", "sha256": ""}]
    job.write_text(json.dumps(artifacts), encoding="utf-8")
    err = _code(lambda: open_dataset(pipeline_path=str(run)))
    assert err.code == "ALIGNMENT_NOT_FOUND"
    assert err.details["candidates"] == ["AlignResult-1_PeakProperties.arf"]


def test_pipeline_folder_falls_back_to_run_input_arf2(tmp_path):
    job = _job(tmp_path, with_folder=False)
    run = _pipeline(tmp_path, job)
    write_msdial_folder(run / "input", **_spec(), stem="AlignResult-1")
    opened = open_dataset(pipeline_path=str(run))
    assert opened.selection.arf2 == run / "input" / "AlignResult-1.arf2"
    assert opened.dataset.source_format == "msdial_folder"


# ---------- ライブラリの記憶 ----------

def test_library_record_per_entry(tmp_path):
    folder = _folder(tmp_path / "gui")
    write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record()], ["k"])])
    assert open_dataset(directory=str(folder)).dataset.library_record["status"] == "known"

    mztab = write_matching_mztab(tmp_path / "direct" / "x.mzTab", _spec())
    assert open_dataset(mztab_path=str(mztab)).dataset.library_record["status"] == "unknown"

    msp = tmp_path / "lab.msp"
    msp.write_text("NAME: x\n", encoding="ascii")
    method = tmp_path / "method.txt"
    method.write_text(f"Msp file path: {msp}\n", encoding="ascii")
    job = _job(tmp_path / "j", method_file=str(method))
    for source in ("auto", "mztab"):
        rec = open_dataset(job_path=str(job), source=source).dataset.library_record
        assert rec["status"] == "known" and rec["origin"] == "effective_method", source

    run_lib = tmp_path / "run_lib.msp"
    run_lib.write_text("NAME: y\n", encoding="ascii")
    run = _pipeline(tmp_path / "p", _job(tmp_path / "p"))
    record_path = run / "pipeline-run.json"
    data = json.loads(record_path.read_text(encoding="utf-8"))
    import hashlib
    data["inputs"] = {"library_mode": "msp_only", "msp": {
        "path": str(run_lib), "file": run_lib.name, "source": "config_file",
        "sha256": hashlib.sha256(run_lib.read_bytes()).hexdigest()}}
    record_path.write_text(json.dumps(data), encoding="utf-8")
    for source in ("auto", "mztab"):
        rec = open_dataset(pipeline_path=str(run), source=source).dataset.library_record
        assert rec["status"] == "known" and rec["origin"] == "run_inputs", source


# ---------- ツール ----------

def test_dataset_load_tool_accepts_directory(tmp_path):
    from metabolomix.tools.handle_tools import load_dataset

    folder = _folder(tmp_path / "gui")
    text = load_dataset(directory=str(folder))
    assert "msdial_folder" in text
    assert f"{_STEM}_PeakProperties.arf" in text
    assert session_state.session.dataset.source_format == "msdial_folder"
    assert str(folder) not in text  # 置き場所は出さない（ファイル名だけ）


def test_dataset_load_tool_reports_library_without_path(tmp_path):
    from metabolomix.tools.handle_tools import load_dataset

    folder = _folder(tmp_path / "gui")
    write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record()], ["k"])])
    text = load_dataset(directory=str(folder))
    assert "Proj_Loaded.msp2.dbs" in text and str(folder) not in text


def test_dataset_load_tool_job_error_names_mztab_source(tmp_path):
    from metabolomix.tools.handle_tools import load_dataset

    job = _job(tmp_path, with_folder=False)
    err = json.loads(load_dataset(job_path=str(job)))["error"]
    assert err["code"] == "ALIGNMENT_NOT_FOUND" and 'source="mztab"' in err["message"]
    assert session_state.session.dataset is None


def test_dataset_handle_open_directory_uses_selection_files(tmp_path):
    from metabolomix.tools import handle_tools

    folder = _folder(tmp_path / "gui")
    (folder / f"{_STEM}.dcl").write_bytes(b"")
    write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record()], ["k"])])
    out = json.loads(handle_tools.load_dataset(directory=str(folder), library="none"))
    assert out["status"] == "success", out
    assert out["adapter"] == "msdial_folder"
    assert out["arf"] == f"{_STEM}_PeakProperties.arf" and out["not_chosen"] == []
    assert out["arf2"] == f"{_STEM}.arf2" and out["dcl"] == f"{_STEM}.dcl"
    assert out["library_record"]["file"] == "Proj_Loaded.msp2.dbs"
    assert "path" not in out["library_record"]
    handle = handles.registry.get(out["dataset_handle"])
    assert handle.arf2_path == str(folder / f"{_STEM}.arf2")
    assert handle.source["adapter"] == "msdial_folder"
    assert str(folder) not in json.dumps(out, ensure_ascii=False)


# ---------- 決めたこと（R15）を縛る ----------

def test_dataset_load_reuses_the_base_and_gives_a_fresh_dataset(tmp_path, monkeypatch):
    """`load_dataset` も基を使い回すが、セッションに置くのは毎回別の実体（別の dataset_id）。"""
    from metabolomix.dataset import msdial_adapter
    from metabolomix.tools.handle_tools import load_dataset

    folder = _folder(tmp_path / "gui")
    calls = []
    real = msdial_adapter.load_msdial_folder

    def spy(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(msdial_adapter, "load_msdial_folder", spy)
    load_dataset(directory=str(folder))
    first = session_state.session.dataset
    first.feature_matrix[0, 0] = -1.0
    load_dataset(directory=str(folder))
    second = session_state.session.dataset
    assert len(calls) == 1
    assert second is not first and second.dataset_id != first.dataset_id
    assert second.feature_matrix[0, 0] == 100.0


@pytest.mark.parametrize("entry", ["mztab_path", "job_path"])
def test_mztab_adapter_rejects_a_non_default_measure_in_both_tools(tmp_path, entry):
    """mzTab-M の定量種別はファイル（ジョブの宣言）に従う。measure は解析フォルダだけに効く。"""
    from metabolomix.tools import handle_tools
    from metabolomix.tools.handle_tools import load_dataset

    job = _job(tmp_path)
    path = str(job.parent / "msdial" / "Height_x.mzTab") if entry == "mztab_path" else str(job)
    kwargs = {entry: path, "measure": "peak_area_above_zero"}
    if entry == "job_path":
        kwargs["source"] = "mztab"
    assert json.loads(load_dataset(**kwargs))["error"]["code"] == "DATASET_BAD_REQUEST"
    out = json.loads(handle_tools.load_dataset(**kwargs, library="none"))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    # 解析フォルダのアダプタでは面積を選べる。
    folder = open_dataset(job_path=str(job), measure="peak_area_above_zero")
    assert folder.dataset.quantification_measure == "peak_area_above_zero"


def test_both_tools_reject_a_non_bool_allow_incomplete(tmp_path):
    """"false" という文字列を真偽として読まない（以前のハンドルは bool() で True にしていた）。"""
    from metabolomix.tools import handle_tools
    from metabolomix.tools.handle_tools import load_dataset

    folder = _folder(tmp_path / "gui")
    load = json.loads(load_dataset(directory=str(folder), allow_incomplete="false"))
    assert load["error"]["code"] == "DATASET_BAD_REQUEST"
    assert session_state.session.dataset is None
    out = json.loads(handle_tools.load_dataset(directory=str(folder), library="none",
                                                      allow_incomplete="false"))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert handles.registry.list() == []


# ---------- 修正 1 回目（レビューの再現を回帰テストにしたもの） ----------

def _other_spec() -> dict:
    """`_spec()` と同じ試料・同じスポット数で、m/z が違う別のアラインメント。"""
    spec = _spec()
    for spot, mz in zip(spec["spots"], (300.1, 400.2)):
        spot["mz"] = mz
    return spec


def test_job_folder_picks_the_batch_that_matches_the_mztab(tmp_path):
    """名前の古いバッチが mzTab-M と合わないなら、合う新しいバッチを選ぶ（先頭を決め打ちしない）。"""
    job = _job(tmp_path)
    write_msdial_folder(job.parent / "msdial", **_other_spec(),
                        stem="AlignmentResult_2026_01_01_00_00_00")
    opened = open_dataset(job_path=str(job))
    assert opened.selection.arf.name == f"{_STEM}_PeakProperties.arf"
    assert opened.dataset.feature_metadata["0"]["mz"] == pytest.approx(760.5851)
    assert any("同じアラインメント" in note for note in opened.notes)


def test_job_folder_with_two_matching_batches_is_ambiguous(tmp_path):
    job = _job(tmp_path)
    write_msdial_folder(job.parent / "msdial", **_spec(), stem="AlignmentResult_2026_03_01_00_00_00")
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert err.code == "ALIGNMENT_AMBIGUOUS"
    assert 'source="mztab"' in err.message
    assert sorted(err.details["candidates"]) == [
        f"{_STEM}_PeakProperties.arf", "AlignmentResult_2026_03_01_00_00_00_PeakProperties.arf"]


def test_job_folder_with_no_matching_batch_is_not_found(tmp_path):
    job = _job(tmp_path, with_folder=False)
    write_msdial_folder(job.parent / "msdial", **_other_spec())
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert err.code == "ALIGNMENT_NOT_FOUND"
    assert 'source="mztab"' in err.message
    assert err.details["candidates"] == [f"{_STEM}_PeakProperties.arf"]


def test_pipeline_input_fallback_is_verified_too(tmp_path):
    job = _job(tmp_path, with_folder=False)
    run = _pipeline(tmp_path, job)
    write_msdial_folder(run / "input", **_other_spec(), stem="AlignResult-1")
    err = _code(lambda: open_dataset(pipeline_path=str(run)))
    assert err.code == "ALIGNMENT_NOT_FOUND"
    assert err.details["candidates"] == ["AlignResult-1_PeakProperties.arf"]


def test_mztab_cache_reloads_when_an_arf_appears(tmp_path):
    """mzTab-M の隣に後から `.arf` が置かれたら、検出状態を取り込み直す（古い基を渡さない）。"""
    folder = tmp_path / "direct"
    mztab = write_matching_mztab(folder / "x.mzTab", _spec())
    assert open_dataset(mztab_path=str(mztab)).dataset.detected_mask is None
    write_msdial_folder(folder, **_spec())
    again = open_dataset(mztab_path=str(mztab))
    assert again.reused_base is False
    assert again.dataset.detected_mask is not None


def test_mztab_cache_reloads_when_the_evidence_arf_changes(tmp_path):
    folder = tmp_path / "direct"
    write_msdial_folder(folder, **_spec())
    mztab = write_matching_mztab(folder / "x.mzTab", _spec())
    first = open_dataset(mztab_path=str(mztab))
    assert first.dataset.detected_mask is not None
    changed = _spec()
    changed["spots"][1]["gap_filled"] = [True, True, False, False]
    write_msdial_folder(folder, **changed)
    again = open_dataset(mztab_path=str(mztab))
    assert again.reused_base is False
    assert again.dataset.detected_mask[1].tolist() == [False, False, True, True]


def test_folder_cache_reloads_when_the_dcl_changes(tmp_path):
    folder = _folder(tmp_path / "gui")
    dcl = folder / f"{_STEM}.dcl"
    dcl.write_bytes(b"")
    open_dataset(directory=str(folder))
    assert open_dataset(directory=str(folder)).reused_base is True
    with open(dcl, "ab") as handle:
        handle.write(b"\0")
    assert open_dataset(directory=str(folder)).reused_base is False


def test_job_cache_reloads_when_the_method_changes(tmp_path):
    """ジョブ経由の読み込みは実効メソッド（と、それが宣言する MSP）も読むので、それも数える。"""
    method = tmp_path / "method.txt"
    method.write_text("Msp file path: \n", encoding="ascii")
    job = _job(tmp_path / "j", method_file=str(method))
    for source in ("auto", "mztab"):
        open_dataset(job_path=str(job), source=source)
        assert open_dataset(job_path=str(job), source=source).reused_base is True, source
    with open(method, "a", encoding="ascii") as handle:
        handle.write("Lbm file path: \n")
    for source in ("auto", "mztab"):
        assert open_dataset(job_path=str(job), source=source).reused_base is False, source


def test_pipeline_cache_reloads_when_the_run_record_changes(tmp_path):
    job = _job(tmp_path)
    run = _pipeline(tmp_path, job)
    open_dataset(pipeline_path=str(run))
    assert open_dataset(pipeline_path=str(run)).reused_base is True
    with open(run / "pipeline-run.json", "a", encoding="utf-8") as handle:
        handle.write("\n")
    assert open_dataset(pipeline_path=str(run)).reused_base is False


def test_folder_polarity_cross_checks_arf2_ion_mode(tmp_path):
    """解析フォルダをジョブ経由で読んでも、宣言極性を `.arf2` の ion_mode の多数で裏取りする。"""
    from metabolomix.mztab.polarity import dataset_polarity

    ds = open_dataset(job_path=str(_job(tmp_path))).dataset  # 宣言は positive
    ds.ion_mode_majority = "negative"
    polarity, origin, notes = dataset_polarity(ds)
    assert (polarity, origin) == ("positive", "job_declared")
    assert any("食い違う" in note and "positive" in note and "negative" in note for note in notes)


def test_job_folder_carries_the_job_artifacts(tmp_path):
    """ジョブ経由で読んだ解析フォルダにも、ジョブが記録した付随アーティファクトを載せる。"""
    from metabolomix.tools.handle_tools import load_dataset

    job = _job(tmp_path)
    data = json.loads(job.read_text(encoding="utf-8"))
    data["artifacts"] = [{"path": "msdial/msdial.log", "role": "console_log",
                          "format": "text", "sha256": ""}]
    job.write_text(json.dumps(data), encoding="utf-8")
    ds = open_dataset(job_path=str(job)).dataset
    assert ds.artifact_paths["console_log"] == [str(job.parent / "msdial" / "msdial.log")]
    assert "console_log" in load_dataset(job_path=str(job))


@pytest.mark.parametrize("entry", ["job_path", "pipeline_path"])
def test_dataset_handle_open_defaults_to_the_folder_for_job_and_pipeline(tmp_path, entry):
    from metabolomix.tools import handle_tools

    job = _job(tmp_path)
    path = job if entry == "job_path" else _pipeline(tmp_path, job)
    out = json.loads(handle_tools.load_dataset(**{entry: str(path)}, library="none"))
    assert out["status"] == "success", out
    assert out["adapter"] == "msdial_folder"
    assert out["arf"] == f"{_STEM}_PeakProperties.arf"


# ---------- 修正 2 回目 ----------

def test_job_folder_cache_reloads_when_the_mztab_changes(tmp_path):
    """解析フォルダの経路も正準 mzTab-M を読んで `.arf` を選ぶので、mzTab-M が変われば読み直す。"""
    job = _job(tmp_path)
    open_dataset(job_path=str(job))
    assert open_dataset(job_path=str(job)).reused_base is True
    write_matching_mztab(job.parent / "msdial" / "Height_x.mzTab", _other_spec())
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert err.code == "ALIGNMENT_NOT_FOUND"


def test_ambiguous_hint_names_directory_only_with_its_caveat(tmp_path):
    """`directory` で開くとジョブの門を通らない。完了したジョブではその注意つきで案内する。"""
    job = _job(tmp_path)
    write_msdial_folder(job.parent / "msdial", **_spec(), stem="AlignmentResult_2026_03_01_00_00_00")
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert err.code == "ALIGNMENT_AMBIGUOUS"
    assert "directory" in err.message and 'source="mztab"' in err.message
    assert "裏取り" in err.message and "探索専用" in err.message


def test_ambiguous_hint_for_an_incomplete_job_names_only_mztab(tmp_path):
    """未完了のジョブで `directory` を案内すると、探索専用の印が外れた読み方を勧めてしまう。"""
    job = _job(tmp_path, status="partial")
    write_msdial_folder(job.parent / "msdial", **_spec(), stem="AlignmentResult_2026_03_01_00_00_00")
    err = _code(lambda: open_dataset(job_path=str(job), allow_incomplete=True))
    assert err.code == "ALIGNMENT_AMBIGUOUS"
    assert 'source="mztab"' in err.message
    assert "directory" not in err.message


def test_mztab_path_with_another_extension_still_tracks_arf_candidates(tmp_path):
    """mzTab-M の読み取りは拡張子を問わない。入力のパスから候補を数えるので、隣に後から
    置かれた `.arf` で読み直す（拡張子で mzTab-M を探すと、ここで候補を見失う）。"""
    folder = tmp_path / "direct"
    mztab = write_matching_mztab(folder / "x.mztab.tsv", _spec())
    first = open_dataset(mztab_path=str(mztab))
    assert first.dataset.detected_mask is None
    write_msdial_folder(folder, **_spec())
    again = open_dataset(mztab_path=str(mztab))
    assert again.reused_base is False
    assert again.dataset.detected_mask is not None


def test_broken_mztab_on_the_folder_route_explains_the_way_out(tmp_path):
    """ジョブの mzTab-M が壊れていると、解析フォルダの経路はバッチを選べない。code は変えずに、
    ジョブの検査なしで解析フォルダを開く道（directory）を案内する。"""
    job = _job(tmp_path)
    (job.parent / "msdial" / "Height_x.mzTab").write_text(
        "MTD\tmzTab-version\t2.0.0-M\n", encoding="utf-8")
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert err.code == "MZTAB_STRUCTURE_INVALID"
    assert "directory" in err.message and "msdial" in err.message


def test_file_changed_during_load_is_not_cached(tmp_path, monkeypatch):
    """読み込みの途中でファイルが変わったら、基をキャッシュしない（次は読み直す）。"""
    from metabolomix.dataset import opening

    folder = _folder(tmp_path / "gui")
    key = f"directory:msdial_folder:peak_height:False:{folder.resolve()}"
    real_open = opening._open

    def changing_open(*args, **kwargs):
        base = real_open(*args, **kwargs)
        # 中身は壊さず、更新時刻だけ進める（読み込み後の読み直しが成り立つように）。
        target = next(folder.glob("*.arf2"))
        stat = target.stat()
        os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns + 5_000_000_000))
        return base

    monkeypatch.setattr(opening, "_open", changing_open)
    first = open_dataset(directory=str(folder))
    assert handles.registry.cached_base(key) is None
    assert any("読み込みの間に入力が変わった" in note for note in first.notes)
    monkeypatch.setattr(opening, "_open", real_open)
    assert open_dataset(directory=str(folder)).reused_base is False


def _bump_mtime_after_load(monkeypatch, target_glob_root: Path, pattern: str):
    """本物の読み込みの後に、`pattern` のファイルの更新時刻だけを進める偽の `_open`。"""
    from metabolomix.dataset import opening

    real_open = opening._open

    def changing_open(*args, **kwargs):
        base = real_open(*args, **kwargs)
        target = next(target_glob_root.rglob(pattern))
        stat = target.stat()
        os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns + 5_000_000_000))
        return base

    monkeypatch.setattr(opening, "_open", changing_open)
    return real_open


def test_job_route_change_during_load_is_not_cached(tmp_path, monkeypatch):
    from metabolomix.dataset import opening

    job = _job(tmp_path)
    real_open = _bump_mtime_after_load(monkeypatch, job.parent / "msdial", "*.arf2")
    first = open_dataset(job_path=str(job))
    assert any("読み込みの間に入力が変わった" in note for note in first.notes)
    monkeypatch.setattr(opening, "_open", real_open)
    assert open_dataset(job_path=str(job)).reused_base is False


def test_non_canonical_path_change_during_load_is_not_cached(tmp_path, monkeypatch):
    from metabolomix.dataset import opening

    folder = _folder(tmp_path / "gui")
    spelled = str(folder).swapcase()  # 大文字小文字違い（Windows ではパスとして同じ）
    if not Path(spelled).exists():
        pytest.skip("大文字小文字を区別するファイルシステム")
    real_open = _bump_mtime_after_load(monkeypatch, folder, "*.arf2")
    first = open_dataset(directory=spelled)
    assert any("読み込みの間に入力が変わった" in note for note in first.notes)
    monkeypatch.setattr(opening, "_open", real_open)
    assert open_dataset(directory=spelled).reused_base is False


def test_pipeline_route_change_during_load_is_not_cached(tmp_path, monkeypatch):
    from metabolomix.dataset import opening

    job = _job(tmp_path)
    run = _pipeline(tmp_path, job)
    real_open = _bump_mtime_after_load(monkeypatch, job.parent / "msdial", "*.arf2")
    first = open_dataset(pipeline_path=str(run))
    assert any("読み込みの間に入力が変わった" in note for note in first.notes)
    monkeypatch.setattr(opening, "_open", real_open)
    assert open_dataset(pipeline_path=str(run)).reused_base is False


# ---------- `.arf` を 1 回だけ読む（B-2a Task 4） ----------

def test_job_folder_reads_arf_once(tmp_path, monkeypatch):
    """job の解析フォルダの経路は、mzTab-M との照合に読んだ `.arf` をそのまま使う（読み直さない）。"""
    from metabolomix.arf import reader as arf_reader

    job = _job(tmp_path)
    calls: list[int] = []
    real = arf_reader.iter_groups

    def counting(file_like):
        calls.append(1)
        return real(file_like)

    monkeypatch.setattr(arf_reader, "iter_groups", counting)
    opened = open_dataset(job_path=str(job))
    assert opened.dataset.source_format == "msdial_folder"
    assert len(opened.dataset.feature_ids) == 2
    assert len(calls) == 1


def test_pipeline_falls_back_to_older_input_batch(tmp_path):
    """run の input/ に複数のバッチがあるとき、最新が合わなくても古いバッチで合えば読む。"""
    job = _job(tmp_path, with_folder=False)
    run = _pipeline(tmp_path, job)
    write_msdial_folder(run / "input", **_other_spec(), stem="AlignResult-2")   # 新しい・合わない
    write_msdial_folder(run / "input", **_spec(), stem="AlignResult-1")         # 古い・合う
    opened = open_dataset(pipeline_path=str(run), source="folder")
    assert opened.selection.arf.name == "AlignResult-1_PeakProperties.arf"
    assert opened.dataset.feature_metadata["0"]["mz"] == pytest.approx(760.5851)


def test_pipeline_arf2_candidates_newest_first(tmp_path):
    from metabolomix.dataset import opening

    job = _job(tmp_path, with_folder=False)
    run = _pipeline(tmp_path, job)
    for stem in ("AlignResult-1", "AlignResult-3", "AlignResult-2"):
        write_msdial_folder(run / "input", **_spec(), stem=stem)
    names = [p.name for p in opening.pipeline_arf2_candidates(run)]
    assert names == ["AlignResult-3.arf2", "AlignResult-2.arf2", "AlignResult-1.arf2"]
    assert opening.pipeline_arf2(run).name == "AlignResult-3.arf2"


def test_broken_mztab_hint_names_reason(tmp_path):
    job = _job(tmp_path)
    (job.parent / "msdial" / "Height_x.mzTab").write_text(
        "MTD\tmzTab-version\t2.0.0-M\n", encoding="utf-8")
    err = _code(lambda: open_dataset(job_path=str(job)))
    assert "読めなかった理由" in err.message and err.code in err.message
    assert 'source="mztab"' in err.message and "でも読めない" in err.message


def test_reports_written_into_the_folder_do_not_invalidate_the_base(tmp_path, monkeypatch):
    """解析フォルダの中に記録（reports/）ができても基を使い回す（フォルダの更新時刻は指紋に入れない）。
    読み込みの最初の概観は記録を reports/stats/ に書くので、これが効かないと 2 回目の読み込みが毎回読み直す。"""
    from metabolomix.dataset import msdial_adapter
    from metabolomix.tools.handle_tools import load_dataset

    folder = _folder(tmp_path / "gui")
    calls = []
    real = msdial_adapter.load_msdial_folder

    def spy(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(msdial_adapter, "load_msdial_folder", spy)
    load_dataset(directory=str(folder), overview=False)
    (folder / "reports" / "stats").mkdir(parents=True)
    (folder / "reports" / "stats" / "x.json").write_text("{}", encoding="utf-8")
    load_dataset(directory=str(folder), overview=False)
    assert len(calls) == 1


def test_folder_listing_counts_mdproject(tmp_path):
    """`.mdproject` は `.mddata` が無いときの群の宣言元。後から置かれたら基を使い回さない（指紋に入る）。"""
    from metabolomix.dataset.opening import _folder_listing

    (tmp_path / "study.mdproject").write_bytes(b"x")
    assert str(tmp_path / "study.mdproject") in _folder_listing(tmp_path)
