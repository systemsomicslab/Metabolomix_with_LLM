"""解析フォルダのアダプタと mzTab-M のアダプタで、同じ job に対するツールの出力がそろうこと
（B-1 最終レビュー R32〜R35）。

アダプタ単体の等価性（`compare_states`）ではなく、ツールが返すもの——検出の要約・キュレーションの
適用・試料のバッチと注入順・実験情報シートの受け付け——を縛る。fixture は tmp に作る。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from metabolomix.core import handles, session_state
from metabolomix.curation import flags as curation_flags
from tests.curation_fixtures import match_result
from tests.msdial_folder_fixtures import analysis_file_row, write_mddata, write_msdial_folder
from tests.test_dataset_opening import _job, _pipeline
from tests.v1_chain import v1_preprocess


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _inchikey(k: int) -> str:
    return "ABCDEFGHIJKLM" + chr(65 + k) + "-OPQRSTUVWX-N"


def _spec() -> dict:
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 3 else "B"}
               for i in range(6)]
    spots = [{"name": f"PC 3{k}:1", "mz": 700.0 + 10 * k,
              "heights": [100.0 + k, 110.0, 120.0, 300.0 + k, 320.0, 340.0],
              "gap_filled": [False] * 5 + [k == 2],
              "msms": [True] * 6,
              "inchikey": _inchikey(k),
              "matches": [match_result({0: f"PC 3{k}:1", 1: _inchikey(k)})]} for k in range(4)]
    return {"spots": spots, "samples": samples}


def _flag_wrong(folder: Path, spot_id: int = 1) -> None:
    arf2 = next(folder.glob("*.arf2"))
    curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
        [{"spot_id": spot_id, "flag": "wrong"}], alignment=curation_flags.alignment_key(arf2),
        review_id="r1", source="test")


def _add_mtd(mztab: Path, batches: list, orders: list) -> None:
    """正準 mzTab-M の MTD に、試料ごとのバッチ（MS:4000088）と注入順（MS:4000089）を足す。"""
    lines = mztab.read_text(encoding="utf-8").splitlines()
    extra = []
    for k, (batch, order) in enumerate(zip(batches, orders), start=1):
        extra.append(f"MTD\tassay[{k}]-custom[1]\t[MS,MS:4000088,batch label,{batch}]")
        extra.append(f"MTD\tassay[{k}]-custom[2]\t[MS,MS:4000089,injection sequence label,{order}]")
    index = max(i for i, line in enumerate(lines) if line.startswith("MTD"))
    lines[index + 1:index + 1] = extra
    mztab.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _load(**kwargs) -> str:
    from metabolomix.tools import handle_tools, mztab_tools

    return handle_tools.load_dataset(**kwargs)


def _sheet(tmp_path: Path, ds, groups: list[str], role: str = "sample") -> Path:
    from metabolomix.analysis import sample_manifest

    paths = [sample_manifest._raw_path_for_sample(ds, i) for i in range(len(ds.sample_names))]
    lines = ["# schema = sample-manifest.v1",
             "sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude"]
    for i, (name, path) in enumerate(zip(ds.sample_names, paths)):
        lines.append(f"{name}\t{Path(path).name}\t{role}\t{groups[i]}\tB9\t{i + 1}\t\ttrue")
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return sheet


GROUPS = ["g1", "g1", "g1", "g2", "g2", "g2"]


# ---------- R32 / 重要 A: 検出状態とキュレーション ----------

def test_folder_dataset_applies_curation_and_reports_detection(tmp_path):
    from metabolomix.tools import dataset_analysis_tools as dat
    from metabolomix.tools import handle_tools, mztab_tools

    folder = tmp_path / "msdial"
    write_msdial_folder(folder, **_spec())
    _flag_wrong(folder)
    _load(directory=str(folder))
    status = json.loads(mztab_tools.dataset_status())
    assert status["detection"]["source"] == "arf"
    assert status["detection"]["n_cells"] == 24 and status["detection"]["n_detected"] == 23
    assert status["inchikey_coverage"]["total_features"] == 4

    v1_preprocess()
    dat.dataset_differential(group_a=["s0", "s1", "s2"], group_b=["s3", "s4", "s5"])
    out = tmp_path / "out.tsv"
    payload = json.loads(dat.dataset_export_differential(str(out), apply_curation=True))
    assert payload["curation"]["state"] == "applied"
    assert payload["curation"]["wrong_excluded"] == 1


# ---------- R33 / 重要 B: 実験情報シート ----------

@pytest.mark.parametrize("entry", ["directory", "job_default", "pipeline"])
def test_set_sample_metadata_accepts_folder_datasets(tmp_path, entry):
    from metabolomix.tools import dataset_analysis_tools as dat

    job = _job(tmp_path, spec=_spec())
    kwargs = {
        "directory": {"directory": str(job.parent / "msdial")},
        "job_default": {"job_path": str(job)},
        "pipeline": {"pipeline_path": str(_pipeline(tmp_path, job))},
    }[entry]
    _load(**kwargs)
    ds = session_state.session.dataset
    assert ds.source_format == "msdial_folder"
    sheet = _sheet(tmp_path, ds, GROUPS)
    result = json.loads(dat.dataset_set_sample_metadata(str(sheet)))
    assert result["status"] == "success", result
    rows = session_state.session.dataset.sample_metadata_rows
    assert [r["group"] for r in rows] == GROUPS
    assert rows[0]["batch"] == "B9"


def test_job_raw_paths_win_over_mddata_file_paths_and_role_is_applied(tmp_path):
    """R25 の優先順: job の mzTab-M の場所が `.mddata` の file_path に勝つ。"""
    from metabolomix.tools import dataset_analysis_tools as dat

    job = _job(tmp_path, spec=_spec())
    mztab = job.parent / "msdial" / "Height_x.mzTab"
    mztab.write_text(mztab.read_text(encoding="utf-8").replace(
        "file:///C:/data/", "file:///D:/elsewhere/"), encoding="utf-8")
    _load(job_path=str(job))
    ds = session_state.session.dataset
    assert ds.msdial_samples[0]["file_path"].replace("\\", "/") == "C:/data/s0.wiff"
    assert ds.sample_raw_paths[0].replace("\\", "/").lower().endswith("d:/elsewhere/s0.wiff")
    # シートは job の場所（D:/elsewhere/）で結ぶ。.mddata の場所では結べない。
    sheet = _sheet(tmp_path, ds, GROUPS, role="qc")
    result = json.loads(dat.dataset_set_sample_metadata(str(sheet)))
    assert result["status"] == "success", result
    rows = session_state.session.dataset.sample_metadata_rows
    assert {r["role"] for r in rows} == {"qc"}
    assert rows[0]["source_file"].replace("\\", "/").lower().endswith("d:/elsewhere/s0.wiff")


def test_set_sample_metadata_without_any_raw_path_says_why(tmp_path):
    from metabolomix.tools import dataset_analysis_tools as dat

    folder = tmp_path / "msdial"
    write_msdial_folder(folder, **_spec(), with_mddata=False)
    _load(directory=str(folder))
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("# schema = sample-manifest.v1\n", encoding="utf-8")
    error = json.loads(dat.dataset_set_sample_metadata(str(sheet)))["error"]
    assert error["code"] == "SAMPLE_MANIFEST_INVALID"
    assert ".mddata" in error["message"] and "console_run" not in error["message"]


# ---------- R34 / 重要 C: MTD のバッチと注入順 ----------

def test_job_folder_carries_mztab_batch_and_injection_order(tmp_path):
    from metabolomix.tools import dataset_analysis_tools as dat

    job = _job(tmp_path, spec=_spec(), with_mddata=False)
    _add_mtd(job.parent / "msdial" / "Height_x.mzTab", [1, 1, 1, 2, 2, 2], [9, 8, 7, 6, 5, 4])
    _load(job_path=str(job))
    ds = session_state.session.dataset
    assert ds.source_format == "msdial_folder"
    assert [ds.sample_meta[n]["batch"] for n in ds.sample_names] == ["1", "1", "1", "2", "2", "2"]
    assert [ds.sample_meta[n]["run_order"] for n in ds.sample_names] == [9, 8, 7, 6, 5, 4]
    assert ds.sample_meta["s0"]["run_order_source"] == "mztab"
    row = ds.sample_metadata_rows[0]
    assert row["provenance"]["injection_order"] == {
        "value": 9, "source": "mztab", "confidence": "unverified"}
    result = json.loads(v1_preprocess())
    assert result["status"] == "success"
    assert "注入順（MTD" not in json.dumps(result, ensure_ascii=False)


def test_mddata_declared_batch_wins_over_mztab_mtd(tmp_path):
    spec = _spec()
    job = _job(tmp_path, spec=spec)
    out = job.parent / "msdial"
    _add_mtd(out / "Height_x.mzTab", [5, 5, 6, 6, 7, 7], [50, 51, 52, 53, 54, 55])
    write_mddata(next(out.glob("*.mddata")), [
        analysis_file_row(s["file_id"], s["name"], s["class_id"], order=100 + s["file_id"],
                          batch=3 + s["file_id"] // 3) for s in spec["samples"]])
    _load(job_path=str(job))
    meta = session_state.session.dataset.sample_meta["s0"]
    assert meta["batch"] == "3" and meta["run_order"] == 100
    assert meta["batch_source"] == "msdial_mddata"


# ---------- R35: 2 つのアダプタでツールの出力がそろう ----------

@pytest.mark.parametrize("entry", ["job", "pipeline"])
def test_folder_and_mztab_adapters_give_the_same_tool_output(tmp_path, entry):
    from metabolomix.analysis import dataset_analysis
    from metabolomix.tools import dataset_analysis_tools as dat
    from metabolomix.tools import handle_tools, mztab_tools

    # mzTab-M のファイル名にバッチ語幹が無くても、m/z の照合で `.arf2`（フラグの置き場所）を見つける。
    job = _job(tmp_path, spec=_spec(), with_mddata=False,
               mztab_name="Height_x.mzTab")
    out = job.parent / "msdial"
    _add_mtd(out / "Height_x.mzTab", [1, 1, 2, 2, 3, 3], [9, 8, 7, 6, 5, 4])
    _flag_wrong(out)
    kwargs = ({"job_path": str(job)} if entry == "job"
              else {"pipeline_path": str(_pipeline(tmp_path, job))})
    seen = {}
    for source in ("folder", "mztab"):
        session_state.session = session_state.AnalysisSession()
        handles.registry = handles.HandleRegistry()
        _load(source=source, **kwargs)
        ds = session_state.session.dataset
        assert ds.source_format == ("msdial_folder" if source == "folder" else "mztab")
        status = json.loads(mztab_tools.dataset_status())
        meta = dataset_analysis.build_dataset_pp_inputs(ds)[-1]
        v1_preprocess()
        dat.dataset_differential(group_a=["s0", "s1", "s2"], group_b=["s3", "s4", "s5"])
        export = json.loads(dat.dataset_export_differential(
            str(tmp_path / f"{source}.tsv"), apply_curation=True))
        applied = json.loads(dat.dataset_set_sample_metadata(str(_sheet(tmp_path, ds, GROUPS))))
        seen[source] = {
            "detection": {k: status["detection"][k]
                          for k in ("available", "source", "n_cells", "n_detected")},
            "curation": (export["curation"]["state"], export["curation"]["wrong_excluded"]),
            "meta": {n: (m["batch"], m["run_order"]) for n, m in meta.items()},
            "applied": applied["status"],
        }
    assert seen["folder"] == seen["mztab"]
    assert seen["folder"]["curation"] == ("applied", 1)
    assert seen["folder"]["meta"]["s0"] == ("1", 9)
    assert seen["folder"]["applied"] == "success"


# ---------- 指摘 5: dataset_status の load_notes と validation ----------

def test_dataset_status_shows_load_notes_and_marks_validation_not_applicable(tmp_path):
    from metabolomix.tools import handle_tools, mztab_tools

    folder = tmp_path / "msdial"
    write_msdial_folder(folder, **_spec(), with_mddata=False)
    _load(directory=str(folder))
    status = json.loads(mztab_tools.dataset_status())
    assert status["validation"]["applicable"] is False and "ok" not in status["validation"]
    assert status["load_notes"]["n"] >= 1
    assert any(".mddata" in line for line in status["load_notes"]["first"])
    assert len(status["load_notes"]["first"]) <= 3


# ---------- R40: 定数のバッチラベルは無視して、ファイル名の日付へ（mzTab 経路と同じ規則） ----------

def _dated_spec() -> dict:
    spec = _spec()
    for sample in spec["samples"]:
        day = "20260101" if sample["file_id"] < 3 else "20260102"
        sample["name"] = f"smp_{day}_{sample['file_id']}"
    return spec


def _differential_caveats(names) -> str:
    from metabolomix.tools import dataset_analysis_tools as dat

    v1_preprocess()
    out = json.loads(dat.dataset_differential(group_a=names[:3], group_b=names[3:]))
    return " ".join(out.get("caveats") or [])


def test_constant_mtd_batch_falls_back_to_filename_dates_on_both_routes(tmp_path):
    from metabolomix.analysis import dataset_analysis

    spec = _dated_spec()
    names = [s["name"] for s in spec["samples"]]
    job = _job(tmp_path, spec=spec, with_mddata=False)
    _add_mtd(job.parent / "msdial" / "Height_x.mzTab", [1] * 6, [9, 8, 7, 6, 5, 4])
    seen = {}
    for source in ("folder", "mztab"):
        session_state.session = session_state.AnalysisSession()
        handles.registry = handles.HandleRegistry()
        _load(job_path=str(job), source=source)
        ds = session_state.session.dataset
        meta = dataset_analysis.build_dataset_pp_inputs(ds)[-1]
        if source == "folder":
            # 行の来歴: バッチは推定（ファイル名の日付）、注入順は mzTab-M の MTD（未確認）。
            prov = ds.sample_metadata_rows[0]["provenance"]
            assert prov["batch"] == {"value": "20260101", "source": "filename_date",
                                     "confidence": "inferred"}
            assert prov["injection_order"]["source"] == "mztab"
        seen[source] = {
            "meta": {n: (m["batch"], m["batch_source"]) for n, m in meta.items()},
            "caveats": _differential_caveats(names)}
    assert seen["folder"]["meta"] == seen["mztab"]["meta"]
    assert seen["folder"]["meta"][names[0]] == ("20260101", "filename_date")
    for source in ("folder", "mztab"):
        assert "交絡:" in seen[source]["caveats"] and "20260101" in seen[source]["caveats"]
        assert "交絡評価不可" not in seen[source]["caveats"]
    assert seen["folder"]["caveats"] == seen["mztab"]["caveats"]


def test_constant_mddata_batch_also_falls_back_to_filename_dates(tmp_path):
    """ARF 経路はバッチを常にファイル名の日付から推定する（MS-DIAL に専用項目が無い）。
    解析フォルダの経路は、定数の `.mddata` のバッチ（GUI の既定の 1）も採らず同じ日付へ戻る。"""
    spec = _dated_spec()
    names = [s["name"] for s in spec["samples"]]
    folder = tmp_path / "msdial"
    write_msdial_folder(folder, **spec)
    _load(directory=str(folder))
    ds = session_state.session.dataset
    assert ds.sample_meta[names[0]]["batch"] == "20260101"
    assert ds.sample_meta[names[0]]["batch_source"] == "filename_date"
    assert "交絡:" in _differential_caveats(names)


# ---------- R42: pipeline の解析フォルダは run の試料対応表の役割を、前処理の前から見せる ----------

def test_pipeline_folder_shows_manifest_roles_in_dataset_status(tmp_path):
    from metabolomix.tools import handle_tools, mztab_tools
    from tests.test_dataset_opening import _manifest_row

    job = _job(tmp_path, spec=_spec())
    rows = [{**_manifest_row(f"s{i}", "g1" if i < 3 else "g2")} for i in range(6)]
    for row in rows:
        row["role"] = "qc" if row["sample_id"] == "s0" else "sample"
        row["provenance"]["role"]["value"] = row["role"]
    _load(pipeline_path=str(_pipeline(tmp_path, job, rows=rows)))
    ds = session_state.session.dataset
    assert ds.source_format == "msdial_folder" and ds.roles["s0"] == "qc"
    assert ds.sample_meta["s0"]["role"] == "qc"
    samples = json.loads(mztab_tools.dataset_status())["samples"]
    assert "s0	qc" in samples and "s1	sample" in samples
