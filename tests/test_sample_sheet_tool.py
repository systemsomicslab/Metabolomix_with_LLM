"""`pipeline_sample_sheet`（ISA / mwTab → sample-manifest.v2）と、`list_data_files` の取得時刻の節。

fixture はテスト自身が tmp に作る（startTimeStamp 付きの合成 mzML）。
"""
from __future__ import annotations

import json
from datetime import timedelta

from metabolomix.analysis import sample_manifest
from metabolomix.pipeline.inputs import list_primary_names
from metabolomix.tools.dataset import list_data_files
from metabolomix.tools.sample_sheet_tools import pipeline_sample_sheet
from tests.test_rawmeta import T0, _isa, write_mzml


def _call(**kw) -> dict:
    return json.loads(pipeline_sample_sheet(**kw))


def test_writes_a_manifest_the_pipeline_reader_accepts(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype", "Diet"],
                assay_table=str(a_table), qc_selector=["Sample Name", "^qc"])
    assert out["status"] == "success", out
    target = raw / "sample-manifest.tsv"
    assert out["path"] == str(target) and target.is_file()
    assert out["injection_order_source"] == "raw_acquisition_time" and out["gaps_over_40min"] == 0
    assert out["groups"] == {"KO_chow": 1, "QC": 1, "WT_chow": 1, "WT_high_fat": 1, "unassigned": 1}
    assert out["columns"][-2:] == ["Genotype", "Diet"]
    rows = sample_manifest.parse_manifest(target, source_root=raw,
                                          expected_sources=list_primary_names(raw, {}))
    by_file = {r["source_file"].replace("\\", "/").rsplit("/", 1)[-1]: r for r in rows}
    assert by_file["r03.mzML"]["injection_order"] == 1 and by_file["qc1.mzML"]["role"] == "qc"
    assert by_file["r01.mzML"]["group"] == "WT_chow" and by_file["r01.mzML"]["Genotype"] == "WT"
    assert any("batch" in note for note in out["notes"])


def test_existing_sheet_is_not_overwritten_unless_asked(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    kw = dict(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"],
              assay_table=str(a_table))
    assert _call(**kw)["status"] == "success"
    again = _call(**kw)
    assert again["error"]["code"] == "SAMPLE_SHEET_EXISTS"
    assert _call(**kw, overwrite=True)["status"] == "success"
    other = _call(**kw, output_name="other.tsv")
    assert other["status"] == "success" and (raw / "other.tsv").is_file()
    bad = _call(**kw, output_name="../x.tsv")
    assert bad["error"]["code"] == "SAMPLE_TABLE_INVALID"


def test_unreadable_raw_leaves_order_empty_and_says_so(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    (raw / "r02.mzML").write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"],
                assay_table=str(a_table))
    assert out["status"] == "success" and out["injection_order_source"] == "none"
    assert out["unreadable"] == ["r02.mzML"]
    assert any("injection_order は空" in note for note in out["notes"])
    rows = sample_manifest.parse_manifest(raw / "sample-manifest.tsv", source_root=raw,
                                          expected_sources=list_primary_names(raw, {}))
    assert all(r["injection_order"] is None for r in rows)


def test_bad_inputs_are_envelopes(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    missing = _call(raw_dir=str(tmp_path / "nope"), sample_table=str(s_table), factor_columns=["Genotype"])
    assert missing["error"]["code"] == "SAMPLE_TABLE_INVALID"
    reserved = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["batch"],
                     assay_table=str(a_table))
    assert reserved["error"]["code"] == "SAMPLE_TABLE_INVALID"
    assert not (raw / "sample-manifest.tsv").exists()


# ---------- list_data_files(all_files=True) の取得時刻 ----------

def test_list_data_files_reports_acquisition_order_and_gaps(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_mzml(raw / "b.mzML", T0)
    write_mzml(raw / "a.mzML", T0 + timedelta(minutes=90))
    (raw / "notes.txt").write_text("x", encoding="utf-8")
    text = list_data_files(directory=str(raw), all_files=True)
    assert "## 取得時刻（raw のメタデータだけを読む。スペクトルは解さない）" in text
    assert "読めた 2 / 2 件" in text and "40 分を超える間隔: 1 箇所" in text
    table = text.split("file\tacquired_utc\torder\tsource\n", 1)[1].splitlines()
    assert table[:2] == [f"a.mzML\t{(T0 + timedelta(minutes=90)).isoformat()}\t2\tmzml",
                         f"b.mzML\t{T0.isoformat()}\t1\tmzml"]
    assert "群と取得順の重なり" not in text                    # シートが無ければ群の判定はしない
    # MS-DIAL の出力だけを見る既定の一覧には取得時刻の節を足さない
    assert "取得時刻" not in list_data_files(directory=str(raw))


def test_list_data_files_reports_group_overlap_from_the_sheet(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    assert _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"],
                 assay_table=str(a_table), qc_selector=["Sample Name", "^qc"])["status"] == "success"
    text = list_data_files(directory=str(raw), all_files=True)
    assert "## 群と取得順の重なり（sample-manifest.tsv の群）" in text
    pairs = text.split("group_a\tgroup_b\toverlap\n", 1)[1].splitlines()
    # 取得順: r03(WT) qc1 r01(WT) r02(KO) r04(unassigned) → WT は 1–3、KO は 4
    assert "KO\tWT\tcomplete" in pairs
    assert "群ごとの取得順の範囲" in text and "QC" not in text.split("群ごとの取得順の範囲", 1)[1]


def test_list_data_files_marks_unreadable_raws(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_mzml(raw / "a.mzML", T0)
    (raw / "b.wiff").write_bytes(b"broken")
    text = list_data_files(directory=str(raw), all_files=True)
    assert "読めた 1 / 2 件" in text
    line = next(x for x in text.splitlines() if x.startswith("b.wiff\t"))
    assert line.split("\t")[2] == "" and line.split("\t")[3].startswith("unreadable")
