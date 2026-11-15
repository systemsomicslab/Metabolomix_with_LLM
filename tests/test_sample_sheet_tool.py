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
from tests.test_rawmeta import T0, _isa, write_agilent_d, write_mzml, write_wiff


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
    # 注入順が取得時刻から来たことは読み直しても残る（出所の名前は raw_acquisition_time。値の出所はシート）
    prov = by_file["r03.mzML"]["provenance"]["injection_order"]
    assert prov["source"] == "user_manifest" and prov["derived_from"] == "raw_acquisition_time"
    assert sample_manifest.order_provenance_label(prov) == "raw_acquisition_time"
    _roles, meta = sample_manifest.pp_inputs_from_rows([r["sample_id"] for r in rows], rows)
    assert {m["run_order_source"] for m in meta.values()} == {"raw_acquisition_time"}


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


def test_unreadable_raw_keeps_the_rest_of_the_order_and_says_which(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    (raw / "r02.mzML").write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"],
                assay_table=str(a_table))
    assert out["status"] == "success" and out["injection_order_source"] == "raw_acquisition_time_partial"
    assert out["unreadable"] == ["r02.mzML"] and out["unreadable_reasons"] == {"開始時刻が無い": 1}
    note = next(n for n in out["notes"] if "取得時刻の読めない raw 1 件" in n)
    assert "ドリフト補正" in note                                   # r02 は生体試料（順が欠ける）
    rows = sample_manifest.parse_manifest(raw / "sample-manifest.tsv", source_root=raw,
                                          expected_sources=list_primary_names(raw, {}))
    order = {r["source_file"].replace("\\", "/").rsplit("/", 1)[-1]: r["injection_order"] for r in rows}
    assert order == {"r03.mzML": 1, "qc1.mzML": 2, "r01.mzML": 3, "r04.mzML": 4, "r02.mzML": None}


def test_no_readable_raw_leaves_every_order_empty(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    for name in ("r01.mzML", "r02.mzML", "r03.mzML", "r04.mzML", "qc1.mzML"):
        (raw / name).write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"],
                assay_table=str(a_table))
    assert out["status"] == "success" and out["injection_order_source"] == "none"
    assert any("injection_order はすべて空" in note for note in out["notes"])
    rows = sample_manifest.parse_manifest(raw / "sample-manifest.tsv", source_root=raw,
                                          expected_sources=list_primary_names(raw, {}))
    assert all(r["injection_order"] is None for r in rows)
    assert all("derived_from" not in r["provenance"]["injection_order"] for r in rows)


def test_agilent_d_folders_get_an_injection_order(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    for k, name in enumerate(["s2.d", "s1.d", "s3.d"]):
        write_agilent_d(raw / name, (T0 + timedelta(minutes=15 * k)).isoformat())
    mwtab = tmp_path / "ST000005.txt"
    mwtab.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000005\n"
                     + "".join(f"SUBJECT_SAMPLE_FACTORS\t-\t{n}\tTreatment:{t}\tRAW_FILE_NAME={n}.d\n"
                               for n, t in (("s1", "ctrl"), ("s2", "drug"), ("s3", "ctrl"))), encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(mwtab), factor_columns=["Treatment"])
    assert out["status"] == "success" and out["injection_order_source"] == "raw_acquisition_time", out
    assert out["unreadable"] == []
    rows = sample_manifest.parse_manifest(raw / "sample-manifest.tsv", source_root=raw,
                                          expected_sources=list_primary_names(raw, {}))
    order = {r["source_file"].replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]: r["injection_order"] for r in rows}
    assert order == {"s2.d": 1, "s1.d": 2, "s3.d": 3}


def test_wiff2_primaries_take_their_order_from_the_same_stem_wiff(tmp_path):
    """SCIEX OS のフォルダ（同じ幹の .wiff・.wiff.scan・.wiff2）で .wiff2 を主ファイルに選んでも、注入順は
    同じ幹の .wiff の取得時刻から付く（mzML への変換を求めない）。list_data_files も .wiff2 を読めたと示す。"""
    raw = tmp_path / "raw"
    raw.mkdir()
    for k, stem in enumerate(["s2", "s1", "s3"]):                  # 1 周期 50 分の取得
        start = T0 + timedelta(minutes=50 * k)
        write_wiff(raw / f"{stem}.wiff", start, start + timedelta(minutes=45))
        (raw / f"{stem}.wiff.scan").write_bytes(b"\x00" * 16)
        (raw / f"{stem}.wiff2").write_bytes(b"\x00" * 64)
    mwtab = tmp_path / "ST000006.txt"
    mwtab.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000006\n"
                     + "".join(f"SUBJECT_SAMPLE_FACTORS\t-\t{n}\tTreatment:{t}\tRAW_FILE_NAME={n}.wiff2\n"
                               for n, t in (("s1", "ctrl"), ("s2", "drug"), ("s3", "drug"))), encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(mwtab), factor_columns=["Treatment"], keep_extension=".wiff2")
    assert out["status"] == "success" and out["injection_order_source"] == "raw_acquisition_time", out
    assert out["unreadable"] == [] and out["unreadable_reasons"] == {}
    names = list_primary_names(raw, {"keep_extension": ".wiff2"})
    rows = sample_manifest.parse_manifest(raw / "sample-manifest.tsv", source_root=raw, expected_sources=names)
    order = {r["source_file"].replace("\\", "/").rsplit("/", 1)[-1]: r["injection_order"] for r in rows}
    assert order == {"s2.wiff2": 1, "s1.wiff2": 2, "s3.wiff2": 3}
    text = list_data_files(directory=str(raw), all_files=True)
    line = next(x for x in text.splitlines() if x.startswith("s1.wiff2\t"))
    assert line.split("\t")[1] == (T0 + timedelta(minutes=50)).isoformat() and line.split("\t")[3] == "wiff_ole_sibling"
    # .wiff と .wiff2 が同じ時刻で並んでも、区切りの間隔は 1 周期（50 分）の 3 倍のまま（間隔 0 で縮めない）
    assert "150 分を超える間隔で区切った取得のまとまり（1 個）" in text


def test_order_missing_only_for_excluded_samples_does_not_claim_drift_is_blocked(tmp_path):
    """取得時刻の読めない raw が統計から外した試料（群 unassigned・include=false）だけなら、ドリフト補正と注入順の
    偏りの判定は止まらない（前処理の _select が先に外す）。notes もそう書き、止まると言わない。"""
    from types import SimpleNamespace

    from metabolomix.analysis.stats import data as stats_data
    from metabolomix.analysis.stats import preprocess as P

    raw = tmp_path / "raw"
    raw.mkdir()
    for k, stem in enumerate(["qc1", "c1", "d1", "qc2", "c2", "d2", "qc3"]):
        write_mzml(raw / f"{stem}.mzML", T0 + timedelta(minutes=10 * k))
    (raw / "x1.mzML").write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")   # 時刻が無く、表に行も無い
    s_table = tmp_path / "s.txt"
    s_table.write_text("Sample Name\tFactor Value[Treatment]\nc1\tCtrl\nc2\tCtrl\nd1\tCase\nd2\tCase\n",
                       encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Treatment"],
                qc_selector=["Sample Name", "^qc"])
    assert out["status"] == "success" and out["unreadable"] == ["x1.mzML"] and out["unassigned_samples"] == ["x1"]
    note = next(n for n in out["notes"] if "取得時刻の読めない raw 1 件" in n)
    assert "その順を埋めるまで行われない" not in note and "止まらない" in note

    names = list_primary_names(raw, {})
    rows = sample_manifest.parse_manifest(raw / "sample-manifest.tsv", source_root=raw, expected_sources=names)

    def drift_reasons(sheet_rows):
        ds = SimpleNamespace(feature_matrix=[[100.0 + k for k in range(len(sheet_rows))]],
                             detected_mask=[[True] * len(sheet_rows)], feature_ids=["f"],
                             sample_names=[r["sample_id"] for r in sheet_rows], sample_metadata_rows=sheet_rows)
        prep = P.preprocess(stats_data.from_dataset(ds), {"filters": {"blank_ratio": None, "qc_rsd": None,
                                                                      "d_ratio": None},
                                                          "normalize": {"method": "none"}, "drift": {"mode": "auto"}})
        step = next(s for s in prep.steps if s["step"] == "drift")
        return step.get("reasons") or []

    assert "no_injection_order" not in drift_reasons(rows)
    # 対照: 同じ試料を統計に入れれば（include=true）、欠けた順がドリフト補正を止める（notes の文言が要る場合）
    included = [dict(r, include=True, group="Ctrl") if r["sample_id"] == "x1" else r for r in rows]
    assert drift_reasons(included) == ["no_injection_order"]


def test_unassigned_renamed_and_ambiguous_rows_are_returned_with_notes(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"],
                assay_table=str(a_table))
    # 因子の値が空の S4 と、表に無い qc1（qc_selector なし）は unassigned・include=false
    assert out["unassigned_samples"] == ["qc1", "S4"] and out["files_not_in_assay"] == ["qc1.mzML"]
    assert any("include=false" in n for n in out["notes"])
    text = (raw / "sample-manifest.tsv").read_text(encoding="utf-8")
    assert "\tunassigned\t\t" in text and all(line.split("\t")[7] == "false" for line in text.splitlines()
                                                if "\tunassigned\t" in line)
    s2 = tmp_path / "s_cols.txt"
    s2.write_text("Sample Name\tbatch\tFactor Value[Genotype]\nS1\t1\tWT\nS2\t1\tKO\nS3\t2\tWT\nS4\t2\tKO\n",
                  encoding="utf-8")
    renamed = _call(raw_dir=str(raw), sample_table=str(s2), factor_columns=["Genotype", "batch"],
                    assay_table=str(a_table), output_name="renamed.tsv")
    assert renamed["status"] == "success" and renamed["renamed_factor_columns"] == {"batch": "factor_batch"}
    assert renamed["columns"][-2:] == ["Genotype", "factor_batch"]
    assert any("'factor_batch'" in n for n in renamed["notes"])


def test_a_sheet_where_nothing_matched_is_an_envelope_and_not_written(tmp_path):
    raw, s_table, _a = _isa(tmp_path)
    wrong = tmp_path / "a_wrong.txt"
    wrong.write_text("Sample Name\tRaw Spectral Data File\nS1\tother_1.raw\nS2\tother_2.raw\n", encoding="utf-8")
    out = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype"], assay_table=str(wrong))
    assert out["error"]["code"] == "SAMPLE_TABLE_INVALID"
    assert out["error"]["details"]["n_files_without_table_row"] == 5
    assert not (raw / "sample-manifest.tsv").exists()


def test_bad_inputs_are_envelopes(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    missing = _call(raw_dir=str(tmp_path / "nope"), sample_table=str(s_table), factor_columns=["Genotype"])
    assert missing["error"]["code"] == "SAMPLE_TABLE_INVALID"
    duplicated = _call(raw_dir=str(raw), sample_table=str(s_table), factor_columns=["Genotype", "genotype"],
                       assay_table=str(a_table))
    assert duplicated["error"]["code"] == "SAMPLE_TABLE_INVALID"
    assert not (raw / "sample-manifest.tsv").exists()


def test_sample_table_that_is_not_utf8_is_an_envelope(tmp_path):
    """latin-1 の ISA の表（µ を含む）は素の例外でなく SAMPLE_TABLE_INVALID の封筒。"""
    raw, _s_table, a_table = _isa(tmp_path)
    latin = tmp_path / "s_latin.txt"
    latin.write_bytes("Sample Name\tFactor Value[Dose]\nS1\t5 \u00b5g\n".encode("latin-1"))
    out = _call(raw_dir=str(raw), sample_table=str(latin), factor_columns=["Dose"], assay_table=str(a_table))
    assert out["error"]["code"] == "SAMPLE_TABLE_INVALID" and "UTF-8" in out["error"]["message"]
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
    # 群 unassigned（因子の値が空の r04）は群として判定しない（組にも範囲にも出ない）
    assert not any("unassigned" in line for line in pairs if "\t" in line)
    assert "unassigned" not in text.split("群ごとの取得順の範囲", 1)[1]
    assert "群 unassigned（群の決まらない試料。include=false）の 1 行は判定に入れない" in text
    # independent の区切りは暦の日付でなく、40 分を超える間隔で区切った取得のまとまり（ここでは 1 個）
    assert "40 分を超える間隔で区切った取得のまとまり（1 個）" in text and "取得日" not in text


def test_list_data_files_marks_unreadable_raws(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_mzml(raw / "a.mzML", T0)
    (raw / "b.wiff").write_bytes(b"broken")
    text = list_data_files(directory=str(raw), all_files=True)
    assert "読めた 1 / 2 件" in text
    line = next(x for x in text.splitlines() if x.startswith("b.wiff\t"))
    assert line.split("\t")[2] == "" and line.split("\t")[3].startswith("unreadable")


def test_list_data_files_counts_raw_and_mzml_twins_as_one_acquisition(tmp_path):
    """raw と変換した mzML が並ぶフォルダでも、1 周期 50 分の取得は 1 つのまとまり（independent）。"""
    from tests.test_rawmeta import write_raw_and_mzml_twins

    raw = tmp_path / "raw"
    raw.mkdir()
    order = ["a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4"]
    write_raw_and_mzml_twins(raw, order, 50)
    (raw / "sample-manifest.tsv").write_text(
        "sample_id\tsource_file\tgroup\trole\n"
        + "".join(f"{n}\t{n}.raw\t{n[0].upper()}\tsample\n" for n in order), encoding="utf-8")
    text = list_data_files(directory=str(raw), all_files=True)
    assert "## 群と取得順の重なり（sample-manifest.tsv の群）: independent" in text
    assert "150 分を超える間隔で区切った取得のまとまり（1 個）" in text
    assert "A\tB\tindependent" in text.split("group_a\tgroup_b\toverlap\n", 1)[1]


def _mwtab(tmp_path, names: list[str], suffix: str):
    table = tmp_path / "ST000007.txt"
    table.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000007\n"
                     + "".join(f"SUBJECT_SAMPLE_FACTORS\t-\t{n}\tTreatment:{'ctrl' if k % 2 else 'drug'}"
                               f"\tRAW_FILE_NAME={n}{suffix}\n" for k, n in enumerate(names)),
                     encoding="utf-8")
    return table


def test_unreadable_reasons_group_files_with_the_same_reason(tmp_path):
    """SCIEX OS のフォルダ（.wiff2 だけで同じ幹の .wiff が無い）でも、理由はファイルによらない 1 種類に
    まとまる（`unreadable` を 20 件に切っても理由が 1 ファイル 1 種類で返り、戻り値が縮まなかった）。"""
    raw = tmp_path / "raw"
    raw.mkdir()
    names = [f"s{k:02d}" for k in range(40)]
    for n in names:
        (raw / f"{n}.wiff2").write_bytes(b"\x00" * 64)
    out = _call(raw_dir=str(raw), sample_table=str(_mwtab(tmp_path, names, ".wiff2")),
                factor_columns=["Treatment"])
    assert out["status"] == "success" and out["injection_order_source"] == "none", out
    assert len(out["unreadable"]) == 20 and out["unreadable_n_more"] == 20
    assert len(out["unreadable_reasons"]) == 1 and list(out["unreadable_reasons"].values()) == [40]
    assert "同じ幹の .wiff" in next(iter(out["unreadable_reasons"]))
    assert "unreadable_reasons_n_more" not in out


def test_unreadable_reasons_are_capped_like_the_lists(tmp_path):
    """ファイルごとに文面の違う理由（壊れた Agilent の AcquiredTime の値）は、件数の多い順に先頭 20 種類だけ返し、
    切った種類の数を `unreadable_reasons_n_more` に置く。"""
    from tests.test_rawmeta import write_agilent_d

    raw = tmp_path / "raw"
    raw.mkdir()
    names = [f"s{k:02d}" for k in range(25)]
    for k, n in enumerate(names):
        # 最後の 3 本は同じ壊れ方（件数 3 の理由が先頭に来る）
        write_agilent_d(raw / f"{n}.d", "broken" if k >= 22 else f"broken{k:02d}")
    out = _call(raw_dir=str(raw), sample_table=str(_mwtab(tmp_path, names, ".d")),
                factor_columns=["Treatment"])
    assert out["status"] == "success", out
    reasons = out["unreadable_reasons"]
    assert len(reasons) == 20 and out["unreadable_reasons_n_more"] == 3     # 23 種類のうち 20 種類
    first, count = next(iter(reasons.items()))
    assert count == 3 and "'broken'" in first
    assert sum(reasons.values()) == 3 + 19
