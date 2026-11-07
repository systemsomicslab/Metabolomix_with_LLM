"""dataset 経路の根拠バンドルが、同じアラインメントの `.dcl` と spectrum_links を書くこと。

結び方は ARF 経路の `alignment_links` と同じ: SMF_ID == AlignmentID == `.dcl` の番号の ID で結び、
m/z と RT は結びを検証するためだけに使う。合わないときは黙って空にせず、manifest の
provenance に状態と理由を書く。
"""
from __future__ import annotations

import csv
import hashlib
import json

import pytest

from metabolomix.analysis import dataset_service
from metabolomix.analysis.evidence_export import export_dataset_evidence
from metabolomix.core import handles, mcp_core, session_state
from metabolomix.mztab import loading
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools
from tests.curation_fixtures import arf2_spot_raw, write_arf2
from tests.dcl_fixture import build_dcl_bytes
from tests.reanalysis_fixtures import POS_FEATURES, write_mztab

HC = ["HC1", "HC2", "HC3"]
RA = ["RA1", "RA2", "RA3"]
EVIDENCE = "same-alignment-stem; SMF_ID == AlignmentID == DCL index; precursor/RT verified"
SPECTRUM = [(60.08, 999.0), (45.03, 120.0)]


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _alignment(run, *, arf2_spots=None, dcl_entries=None):
    """`AlignResult-1.arf2` / `.dcl`。既定は SMF_ID と同じ番号に mzTab と同じ m/z・RT を持つ。

    特徴 1・3 に MS/MS があり、特徴 2 は MS/MS が無い（`no_msms`）。番号 0 は使わない。
    """
    by_id = {int(f["id"]): f for f in POS_FEATURES}
    if arf2_spots is None:
        arf2_spots = [arf2_spot_raw(spot_id=i, name="x", mz=by_id[i]["mz"] if i in by_id else 50.0,
                                    rt=by_id[i]["rt"] if i in by_id else 0.5)
                      for i in range(4)]
    if dcl_entries is None:
        dcl_entries = [{"precursor_mz": 50.0, "rt": 0.5, "spectrum": []}]
        for i in (1, 2, 3):
            dcl_entries.append({"precursor_mz": by_id[i]["mz"], "rt": by_id[i]["rt"],
                                "spectrum": [] if i == 2 else SPECTRUM})
    arf2 = write_arf2(run / "AlignResult-1.arf2", arf2_spots)
    dcl = run / "AlignResult-1.dcl"
    dcl.write_bytes(build_dcl_bytes(dcl_entries))
    return arf2, dcl


def _inputs(tmp_path):
    ds = loading.load_dataset_state(mztab_path=write_mztab(tmp_path / "run" / "AlignResult-1.mzTab"))
    dataset_service.preprocess_dataset(ds, {"normalize": "none"})
    result = dataset_service.compare_dataset(ds, HC, RA, group_a_label="HC", group_b_label="RA")
    return ds, result


def _manifest(out):
    return json.loads((out / "manifest.json").read_text(encoding="utf-8"))


def test_links_are_written_by_id_and_dcl_is_copied(tmp_path):
    ds, result = _inputs(tmp_path)
    arf2, dcl = _alignment(tmp_path / "run")
    out = tmp_path / "b"
    export_dataset_evidence(ds, result, out, polarity="positive", arf2_path=arf2, dcl_path=dcl)
    manifest = _manifest(out)
    source_id = manifest["source_analysis_id"]
    links = _tsv(out / "spectrum_links.tsv")
    assert links == [dict(feature_id=fid, spectrum_file="AlignResult-1.dcl",
                          spectrum_id=f"AlignResult-1:{fid}", source_analysis_id=source_id,
                          mapping_evidence=EVIDENCE) for fid in ("1", "3")]
    assert (out / "AlignResult-1.dcl").read_bytes() == dcl.read_bytes()
    spectra = [f for f in manifest["files"] if f["role"] == "spectra"]
    assert spectra == [dict(role="spectra", path="AlignResult-1.dcl",
                            sha256=hashlib.sha256(dcl.read_bytes()).hexdigest())]
    prov = manifest["provenance"]
    assert prov["spectrum_links_status"] == "verified_alignment"
    assert prov["spectrum_links_reason"] is None
    assert prov["alignment_files"] == {"arf2": "AlignResult-1.arf2", "dcl": "AlignResult-1.dcl"}
    assert prov["msms_link_statuses"] == {"1": {"status": "ok", "dcl_index": 1},
                                          "2": {"status": "no_msms", "dcl_index": 2},
                                          "3": {"status": "ok", "dcl_index": 3}}
    assert prov["msms_link_status_counts"] == {"ok": 2, "no_msms": 1}
    assert prov["id_space_check"] == {"n_checked": 3, "n_mismatch": 0, "examples": []}


def test_missing_dcl_is_recorded_not_silently_empty(tmp_path):
    ds, result = _inputs(tmp_path)
    arf2, _ = _alignment(tmp_path / "run")
    out = tmp_path / "b"
    export_dataset_evidence(ds, result, out, arf2_path=arf2, dcl_path=None)
    prov = _manifest(out)["provenance"]
    assert _tsv(out / "spectrum_links.tsv") == []
    assert not list(out.glob("*.dcl"))
    assert prov["spectrum_links_status"] == "not_provided"
    assert "dcl" in prov["spectrum_links_reason"]
    assert prov["alignment_files"] == {"arf2": "AlignResult-1.arf2", "dcl": None}


def test_no_alignment_files_records_reason(tmp_path):
    ds, result = _inputs(tmp_path)
    out = tmp_path / "b"
    export_dataset_evidence(ds, result, out)
    prov = _manifest(out)["provenance"]
    assert prov["spectrum_links_status"] == "not_provided"
    assert "arf2" in prov["spectrum_links_reason"]


def test_id_space_mismatch_writes_no_links_and_says_why(tmp_path):
    """AlignmentID n の m/z が SMF_ID n の m/z と違う（番号がずれた別のアラインメント）。"""
    ds, result = _inputs(tmp_path)
    shifted = [arf2_spot_raw(spot_id=i, name="x", mz=mz, rt=rt)
               for i, (mz, rt) in enumerate([(104.1070, 1.20), (184.0733, 1.25),
                                             (500.2000, 5.00), (700.0, 7.0)])]
    arf2, dcl = _alignment(tmp_path / "run", arf2_spots=shifted)
    out = tmp_path / "b"
    export_dataset_evidence(ds, result, out, arf2_path=arf2, dcl_path=dcl)
    manifest = _manifest(out)
    prov = manifest["provenance"]
    assert _tsv(out / "spectrum_links.tsv") == []
    assert not [f for f in manifest["files"] if f["role"] == "spectra"]
    assert prov["spectrum_links_status"] == "id_space_mismatch"
    assert "SMF_ID" in prov["spectrum_links_reason"]
    assert prov["id_space_check"]["n_checked"] == 3
    assert prov["id_space_check"]["n_mismatch"] == 3
    assert prov["id_space_check"]["examples"][0]["feature_id"] == "1"


def test_dcl_precursor_mismatch_is_per_feature(tmp_path):
    """`.dcl` の番号 3 の precursor が `.arf2` と合わない → その特徴だけ結ばず、理由を残す。"""
    ds, result = _inputs(tmp_path)
    entries = [{"precursor_mz": 50.0, "rt": 0.5, "spectrum": []},
               {"precursor_mz": 104.1070, "rt": 1.20, "spectrum": SPECTRUM},
               {"precursor_mz": 184.0733, "rt": 1.25, "spectrum": []},
               {"precursor_mz": 640.0, "rt": 5.00, "spectrum": SPECTRUM}]
    arf2, dcl = _alignment(tmp_path / "run", dcl_entries=entries)
    out = tmp_path / "b"
    export_dataset_evidence(ds, result, out, arf2_path=arf2, dcl_path=dcl)
    prov = _manifest(out)["provenance"]
    assert [r["feature_id"] for r in _tsv(out / "spectrum_links.tsv")] == ["1"]
    assert prov["spectrum_links_status"] == "verified_alignment"
    assert prov["msms_link_statuses"]["3"] == {"status": "precursor_mismatch", "dcl_index": 3}
    assert prov["msms_link_status_counts"] == {"ok": 1, "no_msms": 1, "precursor_mismatch": 1}


def test_dcl_shorter_than_ids_is_out_of_range(tmp_path):
    ds, result = _inputs(tmp_path)
    entries = [{"precursor_mz": 50.0, "rt": 0.5, "spectrum": []},
               {"precursor_mz": 104.1070, "rt": 1.20, "spectrum": SPECTRUM}]
    arf2, dcl = _alignment(tmp_path / "run", dcl_entries=entries)
    out = tmp_path / "b"
    export_dataset_evidence(ds, result, out, arf2_path=arf2, dcl_path=dcl)
    statuses = _manifest(out)["provenance"]["msms_link_statuses"]
    assert statuses["2"] == {"status": "out_of_range"}
    assert statuses["3"] == {"status": "out_of_range"}


def _differential(handle=None):
    T.dataset_preprocess(dataset_handle=handle)
    return json.loads(T.dataset_differential(group_a=HC, group_b=RA, dataset_handle=handle))["result_id"]


def test_tool_on_handle_uses_the_handles_alignment_files(tmp_path):
    path = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab")
    _alignment(tmp_path / "run")
    opened = json.loads(handle_tools.load_dataset(mztab_path=str(path), library="none"))
    handle = opened["dataset_handle"]
    result_id = _differential(handle)
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id,
                                                      dataset_handle=handle))
    assert out["status"] == "success", out
    assert [r["feature_id"] for r in _tsv(tmp_path / "b" / "spectrum_links.tsv")] == ["1", "3"]
    assert _manifest(tmp_path / "b")["provenance"]["spectrum_links_status"] == "verified_alignment"


def test_tool_on_session_finds_the_alignment_files(tmp_path, monkeypatch):
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    path = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab")
    _alignment(tmp_path / "run")
    from metabolomix.tools import mztab_tools
    handle_tools.load_dataset(mztab_path=str(path))
    result_id = _differential()
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id))
    assert out["status"] == "success", out
    assert [r["feature_id"] for r in _tsv(tmp_path / "b" / "spectrum_links.tsv")] == ["1", "3"]
