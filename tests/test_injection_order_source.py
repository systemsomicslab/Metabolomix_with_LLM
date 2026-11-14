"""注入順の出所と「MS-DIAL の既定（ファイルの順）かもしれない」印（A6）、シートの空欄の .mddata での補い。

対象: `sample_manifest` の出所の部品と `_merge_row`・シートの宣言の行、交絡の判定（`confounding_report` と、
`dataset_compare` / analysis-spec v2 の比較が組む交絡の判定）、質量精度の注入順の傾向（`mass_accuracy.injection_trend`）。
fixture はテスト自身が作る。
"""
from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from metabolomix.analysis import mass_accuracy, sample_manifest
from metabolomix.analysis.stats import data as stats_data
from metabolomix.analysis.stats import preprocess as P
from metabolomix.analysis.stats.confounding import confounding_report
from metabolomix.core import session_state
from metabolomix.core.atomic_io import DomainError


def test_matches_file_order_rules():
    assert sample_manifest.matches_file_order([2, 1, 3], ["b", "a", "c"])
    assert sample_manifest.matches_file_order([1, 2, 3], ["A1", "a2", "B3"])          # 大文字小文字を無視
    assert not sample_manifest.matches_file_order([1, 3, 2], ["a", "b", "c"])
    assert not sample_manifest.matches_file_order([1, 2], ["a", "b"])                  # 3 試料未満は判定しない
    assert sample_manifest.matches_file_order([1, None, 2, 3], ["a", "z", "b", "c"])   # 値の無い試料は数えない


def test_matches_file_order_rejects_degenerate_orders():
    """すべて同じ・同順位を含む注入順は MS-DIAL の既定（1, 2, …）ではない。同順位をファイル名で割って
    「一致」にしない（drift=auto を保留し、補正できない drift_mode="on" を勧めてしまう）。"""
    assert not sample_manifest.matches_file_order([1, 1, 1], ["c", "a", "b"])
    assert not sample_manifest.matches_file_order([5, 5, 5, 5], ["d", "b", "a", "c"])
    assert not sample_manifest.matches_file_order([1, 1, 2, 3], ["a", "b", "c", "d"])
    assert sample_manifest.matches_file_order([10, 20, 30], ["a", "b", "c"])           # 間が空いても狭義に増えればよい


def test_source_labels():
    label = sample_manifest.injection_order_source_label
    assert [label(s) for s in (None, "user_manifest", "msdial_mddata", "mztab",
                               "mztab_injection_sequence", "raw_acquisition_time")] == [
        "sheet", "sheet", "mddata", "mzTab", "mzTab", "raw_acquisition_time"]


def _msdial_ds(orders, *, batches=("1", "1", "1"), names=("s_c", "s_a", "s_b")):
    samples = [{"name": n, "role": "sample", "group": "G", "batch": b, "injection_order": o,
                "file_path": f"C:/raw/{n}.wiff"} for n, b, o in zip(names, batches, orders)]
    return SimpleNamespace(sample_names=list(names), msdial_samples=samples, sample_assay_ids=[],
                           assay_metadata={}, assay_sources={})


def test_file_order_warning_from_default_rows():
    ds = _msdial_ds([3, 1, 2])                      # s_a=1, s_b=2, s_c=3 → ファイル名の順
    warn = sample_manifest.injection_order_file_order_warning(ds)
    assert warn["code"] == "injection_order_may_be_file_order" and warn["n_samples"] == 3
    assert sample_manifest.injection_order_file_order_warning(_msdial_ds([1, 2, 3])) is None


def _sheet_row(name, *, batch=None, order=None):
    explicit = {"source": "user_manifest", "confidence": "confirmed"}
    empty = {"value": None, "source": None, "confidence": None}
    return {"sample_id": name, "source_file": f"C:/raw/{name}.wiff", "role": "sample", "group": "G",
            "batch": batch, "injection_order": order, "qc_pool": None, "include": True,
            "biological_sample_id": None,
            "provenance": {"sample_id": {"value": name, **explicit},
                           "source_file": {"value": f"C:/raw/{name}.wiff", **explicit},
                           "role": {"value": "sample", **explicit}, "group": {"value": "G", **explicit},
                           "batch": {"value": batch, **explicit} if batch else dict(empty),
                           "injection_order": {"value": order, **explicit} if order else dict(empty),
                           "qc_pool": dict(empty), "include": {"value": True, **explicit},
                           "biological_sample_id": dict(empty)},
            "conflicts": []}


def test_sheet_blanks_are_filled_from_mddata_with_source():
    """A6: シートの空欄の注入順・バッチは .mddata の宣言で埋め、値と出所を行に残す（統計と同じ値）。"""
    ds = _msdial_ds([3, 1, 2], batches=("X", "Y", "Y"))
    merged = sample_manifest._merge_row(ds, 0, "s_c", _sheet_row("s_c"))
    assert merged["injection_order"] == 3 and merged["batch"] == "X"
    assert merged["provenance"]["injection_order"] == {"value": 3, "source": "msdial_mddata",
                                                       "confidence": "declared"}
    assert merged["provenance"]["batch"]["source"] == "msdial_mddata"
    # シートの値は上書きせず、食い違う .mddata の値は conflicts に残す
    kept = sample_manifest._merge_row(ds, 0, "s_c", _sheet_row("s_c", batch="B9", order=7))
    assert kept["injection_order"] == 7 and kept["batch"] == "B9"
    assert {c["field"] for c in kept["conflicts"]} == {"batch", "injection_order"}
    # 定数の .mddata のバッチ（GUI の既定 1）は情報が無いので埋めない
    const = sample_manifest._merge_row(_msdial_ds([3, 1, 2]), 0, "s_c", _sheet_row("s_c"))
    assert const["batch"] is None and const["injection_order"] == 3


def test_default_row_batch_keeps_filename_date_only_as_candidate():
    ds = _msdial_ds([1, 2, 3], names=("smp_20260101_a", "smp_20260101_b", "smp_20260102_c"))
    rows = sample_manifest.resolve_metadata(ds, None)
    assert [r["batch"] for r in rows] == [None, None, None]
    assert rows[2]["provenance"]["batch"]["candidate"] == {"value": "20260102", "source": "filename_date",
                                                           "confidence": "inferred"}


def test_confounding_report_flags_file_order():
    units = [SimpleNamespace(group=g, batch="b1", injection_order=float(o),
                             injection_order_may_be_file_order=True)
             for g, o in (("A", 1), ("A", 2), ("A", 3), ("B", 4), ("B", 5), ("B", 6))]
    rep = confounding_report(units, ["A", "B"])
    assert rep["injection_order"]["may_be_file_order"] is True
    warn = next(w for w in rep["warnings"] if w["code"] == "injection_order_may_be_file_order")
    assert warn["n_units"] == 6
    plain = [SimpleNamespace(group=u.group, batch=u.batch, injection_order=u.injection_order) for u in units]
    rep = confounding_report(plain, ["A", "B"])
    assert "may_be_file_order" not in rep["injection_order"]
    assert not any(w["code"] == "injection_order_may_be_file_order" for w in rep["warnings"])


# ---------- 比較（dataset_compare・analysis-spec v2）の交絡の判定 ----------

def _file_order_ds(tmp_path):
    """試料の表の注入順が .mddata 由来（MS-DIAL の宣言）で、試料名（ファイルの幹）の順と一致する DatasetState。"""
    from tests.test_stats_results import _ds

    ds = _ds()
    rank = {n: i for i, n in enumerate(sorted(ds.sample_names), start=1)}
    for name, row in zip(ds.sample_names, ds.sample_metadata_rows):
        row["injection_order"] = rank[name]
        row["provenance"] = {"injection_order": {"value": rank[name], "source": "msdial_mddata",
                                                 "confidence": "declared"}}
    ds.analysis_dir = str(tmp_path / "analysis")
    return ds


def test_dataset_compare_confounding_carries_the_file_order_warning(tmp_path):
    """A6: 比較の交絡の判定にも「注入順は MS-DIAL の既定（ファイルの順）かもしれない」が届く。"""
    from metabolomix.tools.stats_tools import dataset_compare, dataset_preprocess

    ds = _file_order_ds(tmp_path)
    assert all(s.injection_order_may_be_file_order for s in stats_data.from_dataset(ds).samples)
    session_state.session.dataset = ds
    pre = json.loads(dataset_preprocess(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                                        max_qc_rsd=None, max_d_ratio=None, normalize="none"))
    assert pre["status"] == "success", pre
    out = json.loads(dataset_compare(reference_group="C", test_group="T", matrix_id=pre["matrix_id"]))
    assert out["status"] == "success", out
    conf = [w for w in out["confounding"]["warnings"] if w["code"] == "injection_order_may_be_file_order"]
    assert len(conf) == 1 and conf[0]["n_units"] == 8
    # 前処理の注意（n_samples）は上の warnings に 1 つだけ（交絡の判定の写しを重ねない）
    top = [w for w in out["warnings"] if w["code"] == "injection_order_may_be_file_order"]
    assert len(top) == 1 and "n_samples" in top[0]
    engine = ds.results[out["result_id"]]["engine"]
    assert engine["confounding"]["injection_order"]["may_be_file_order"] is True


def test_spec_v2_compare_confounding_carries_the_file_order_flag(tmp_path):
    from tests.test_stats_results import T_VS_C, _run

    _n, _prepared, result = _run(_file_order_ds(tmp_path), [T_VS_C])
    conf = result["analyses"][0]["result"]["confounding"]
    assert conf["injection_order"]["may_be_file_order"] is True
    assert any(w["code"] == "injection_order_may_be_file_order" for w in conf["warnings"])


# ---------- シートの宣言の行（注入順を何から作ったか） ----------

_V2_HEADER = "sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude\tbiological_sample_id"


def _declared_sheet(tmp_path, declarations: list[str], orders=("2", "1", "3")):
    root = tmp_path / "raw"
    root.mkdir(exist_ok=True)
    names = ["r1.mzML", "r2.mzML", "r3.mzML"]
    for name in names:
        (root / name).write_bytes(b"x")
    body = [f"s{i}\t{n}\tsample\tG\t\t{o}\t\ttrue\t" for i, (n, o) in enumerate(zip(names, orders))]
    sheet = root / "sample-manifest.tsv"
    sheet.write_text("\n".join(["# schema = sample-manifest.v2", *declarations, _V2_HEADER, *body]) + "\n",
                     encoding="utf-8")
    return sheet, root, names


@pytest.mark.parametrize("declared", ["raw_acquisition_time", "filename"])
def test_declared_order_source_is_kept_as_derived_from(tmp_path, declared):
    sheet, root, names = _declared_sheet(tmp_path, [f"# injection_order_source = {declared}"])
    rows = sample_manifest.parse_manifest(sheet, source_root=root, expected_sources=names)
    prov = rows[0]["provenance"]["injection_order"]
    assert prov == {"value": 2, "source": "user_manifest", "confidence": "confirmed", "derived_from": declared}
    assert sample_manifest.order_provenance_label(prov) == declared
    _roles, meta = sample_manifest.pp_inputs_from_rows([r["sample_id"] for r in rows], rows)
    assert {m["run_order_source"] for m in meta.values()} == {declared}
    # 統計エンジンの試料の出所・試料の表の出所の欄も宣言の名前
    ds = SimpleNamespace(feature_matrix=[[1.0, 2.0, 3.0]], detected_mask=None, feature_ids=["f1"],
                         sample_names=[r["sample_id"] for r in rows], sample_metadata_rows=rows)
    assert [s.injection_order_source for s in stats_data.from_dataset(ds).samples] == [declared] * 3
    from metabolomix.tools.mztab_tools import _source_cell
    assert _source_cell(rows[0], "injection_order") == declared


def test_sheet_declaration_and_blank_orders(tmp_path):
    """`sheet` の宣言は宣言しないのと同じ。値の無い欄には derived_from を足さない。"""
    sheet, root, names = _declared_sheet(tmp_path, ["# injection_order_source = sheet"])
    rows = sample_manifest.parse_manifest(sheet, source_root=root, expected_sources=names)
    assert all("derived_from" not in r["provenance"]["injection_order"] for r in rows)
    assert sample_manifest.order_provenance_label(rows[0]["provenance"]["injection_order"]) == "sheet"
    sheet, root, names = _declared_sheet(tmp_path, ["# injection_order_source = raw_acquisition_time"],
                                         orders=("", "", ""))
    rows = sample_manifest.parse_manifest(sheet, source_root=root, expected_sources=names)
    assert all(r["provenance"]["injection_order"] == {"value": None, "source": None, "confidence": None}
               for r in rows)


def test_blank_orders_of_a_partly_ordered_declared_sheet_are_not_filled_from_msdial(tmp_path):
    """取得時刻の読めない raw の空欄（一部だけ順のある宣言つきのシート）は、.mddata の番号で埋めない。

    シートの順位（1, 2）と MS-DIAL の番号（10, 20, 30）は付け方が違う。埋めると r2 が 20 になり、ドリフト補正が
    r2 を r1・r3 より後ろに置いてしまう。空欄は不明（None）のままにし、ドリフト補正は順の欠けで止まる。
    """
    sheet, root, names = _declared_sheet(tmp_path, ["# injection_order_source = raw_acquisition_time"],
                                         orders=("2", "", "1"))
    rows = sample_manifest.parse_manifest(sheet, source_root=root, expected_sources=names)
    assert [sample_manifest.sheet_order_unknown(r) for r in rows] == [False, True, False]
    samples = [{"name": f"r{i}", "role": "sample", "group": "G", "batch": "1", "injection_order": o,
                "file_path": str(root / n)} for i, (n, o) in enumerate(zip(names, (10, 20, 30)))]
    ds = SimpleNamespace(sample_names=[s["name"] for s in samples], msdial_samples=samples, sample_assay_ids=[],
                         assay_metadata={}, assay_sources={}, feature_matrix=[[1.0, 2.0, 3.0]],
                         detected_mask=None, feature_ids=["f1"])
    merged = sample_manifest.resolve_metadata(ds, rows)
    assert [r["injection_order"] for r in merged] == [2, None, 1]
    ds.sample_metadata_rows = merged
    study = stats_data.from_dataset(ds)
    assert [s.injection_order for s in study.samples] == [2.0, None, 1.0]
    assert [s.injection_order_source for s in study.samples] == ["raw_acquisition_time", None,
                                                                 "raw_acquisition_time"]
    _roles, meta = sample_manifest.pp_inputs_from_rows(ds.sample_names, merged)
    assert [m["run_order_source"] for m in meta.values()] == ["raw_acquisition_time", None, "raw_acquisition_time"]
    # 宣言の無いシート（利用者が手で書いた）の空欄は、これまでどおり .mddata で埋める
    (tmp_path / "plain").mkdir()
    plain, root2, names2 = _declared_sheet(tmp_path / "plain", [], orders=("2", "", "1"))
    plain_rows = sample_manifest.parse_manifest(plain, source_root=root2, expected_sources=names2)
    for s, n in zip(samples, names2):
        s["file_path"] = str(root2 / n)
    assert [r["injection_order"] for r in sample_manifest.resolve_metadata(ds, plain_rows)] == [2, 20, 1]


@pytest.mark.parametrize("lines", [["# injection_order_source = mddata"], ["# note = x"], ["# free comment"],
                                   ["# injection_order_source = filename",
                                    "# injection_order_source = filename"]])
def test_unknown_or_duplicated_declarations_are_rejected(tmp_path, lines):
    sheet, root, names = _declared_sheet(tmp_path, lines)
    with pytest.raises(DomainError) as e:
        sample_manifest.parse_manifest(sheet, source_root=root, expected_sources=names)
    assert e.value.code == "SAMPLE_MANIFEST_INVALID"


def test_filename_declared_order_holds_auto_drift():
    """シートの宣言でファイル名の順の注入順は確かでない: drift=auto は補正しない（on なら補正する）。"""
    from dataclasses import replace

    from tests.test_stats_preprocess import _drift_study

    study = _drift_study()
    study.samples = [replace(s, injection_order_source="filename") for s in study.samples]
    cfg = {"filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None},
           "normalize": {"method": "none"}}
    prep = P.preprocess(study, cfg)
    drift = next(s for s in prep.steps if s["step"] == "drift")
    assert drift["status"] == "skipped" and drift["reasons"] == ["injection_order_may_be_file_order"]
    on = P.preprocess(study, {**cfg, "drift": {"mode": "on"}})
    assert next(s for s in on.steps if s["step"] == "drift")["status"] == "applied"


def test_mass_accuracy_rho_carries_the_file_order_warning():
    names = [f"s{i}" for i in range(6)]
    ds = _msdial_ds(list(range(1, 7)), batches=("1",) * 6, names=tuple(names))
    ds.sample_metadata_rows = sample_manifest.resolve_metadata(ds, None)
    ds.per_sample_mz = {"f1": {n: 100.0 + 1e-4 * i for i, n in enumerate(names)}}
    rows = [{"feature_id": "f1", "theoretical_mz": 100.0}]
    trend = mass_accuracy.injection_trend(ds, rows)
    assert trend["injection_order_source"] == "metadata" and trend["spearman_rho"] is not None
    assert trend["warnings"][0]["code"] == "injection_order_may_be_file_order"
    shuffled = _msdial_ds([2, 1, 4, 3, 6, 5], batches=("1",) * 6, names=tuple(names))
    shuffled.sample_metadata_rows = sample_manifest.resolve_metadata(shuffled, None)
    shuffled.per_sample_mz = ds.per_sample_mz
    assert "warnings" not in mass_accuracy.injection_trend(shuffled, rows)


# ---------- シートの宣言でファイル名の順（`# injection_order_source = filename`）の注入順 ----------

def _filename_declared_ds(tmp_path, declaration="filename"):
    """6 試料（A,A,A,B,B,B）・注入順 1..6 を、シートが宣言でファイル名の順と書いた（または宣言なしの）データセット。"""
    root = tmp_path / f"raw_{declaration or 'plain'}"
    root.mkdir()
    names = [f"r{i}.mzML" for i in range(1, 7)]
    for name in names:
        (root / name).write_bytes(b"x")
    body = [f"s{i}\t{n}\tsample\t{'A' if i <= 3 else 'B'}\t\t{i}\t\ttrue\t"
            for i, n in enumerate(names, start=1)]
    head = [f"# injection_order_source = {declaration}"] if declaration else []
    sheet = root / "sample-manifest.tsv"
    sheet.write_text("\n".join(["# schema = sample-manifest.v2", *head, _V2_HEADER, *body]) + "\n",
                     encoding="utf-8")
    rows = sample_manifest.parse_manifest(sheet, source_root=root, expected_sources=names)
    sample_ids = [r["sample_id"] for r in rows]
    return SimpleNamespace(feature_matrix=[[float(10 + i) for i in range(6)]], detected_mask=None,
                           feature_ids=["f1"], sample_names=sample_ids, sample_metadata_rows=rows,
                           per_sample_mz={"f1": {n: 100.0 + 1e-4 * k for k, n in enumerate(sample_ids)}})


def test_filename_declared_order_suppresses_mass_accuracy_rho(tmp_path):
    """A6: シートが宣言でファイル名の順と書いた注入順は、自前のファイル名の順と同じく ρ を出さない
    （試料の表の出所 `filename` と、質量精度・割れの組の判定の出所を揃える）。"""
    from metabolomix.analysis import alignment_qc

    ds = _filename_declared_ds(tmp_path)
    order, source = alignment_qc.sample_injection_order(ds)
    assert source == "filename" and order == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]   # 値はシートのまま
    rows = [{"feature_id": "f1", "theoretical_mz": 100.0}]
    trend = mass_accuracy.injection_trend(ds, rows)
    assert trend["injection_order_source"] == "filename"
    assert trend["spearman_rho"] is None and trend["rho_unavailable_reason"] == "filename_order"
    # 宣言の無いシート（利用者が手で入れた順）は従来どおり manifest で ρ を出す
    plain = _filename_declared_ds(tmp_path, declaration=None)
    assert alignment_qc.sample_injection_order(plain)[1] == "manifest"
    assert mass_accuracy.injection_trend(plain, rows)["spearman_rho"] == pytest.approx(1.0)


def test_filename_declared_order_carries_the_file_order_caveat_in_confounding(tmp_path):
    """A6: 交絡の判定も、宣言でファイル名の順の注入順の偏りに「ファイル名の並びかもしれない」を添える。"""
    from scipy import stats as sps

    ds = _filename_declared_ds(tmp_path)
    samples = stats_data.from_dataset(ds).samples
    assert [s.injection_order_source for s in samples] == ["filename"] * 6
    assert not any(s.injection_order_may_be_file_order for s in samples)    # matrix_id を変えない
    rep = confounding_report(samples, ["A", "B"])
    expected_p = sps.kruskal([1, 2, 3], [4, 5, 6]).pvalue                   # 独立の計算
    assert rep["injection_order"]["status"] == "biased"
    assert rep["injection_order"]["kruskal_p"] == pytest.approx(expected_p)
    assert rep["injection_order"]["may_be_file_order"] is True
    warn = next(w for w in rep["warnings"] if w["code"] == "injection_order_may_be_file_order")
    assert warn["n_units"] == 6 and "injection_order_source = filename" in warn["message"]
    plain = stats_data.from_dataset(_filename_declared_ds(tmp_path, declaration=None)).samples
    rep = confounding_report(plain, ["A", "B"])
    assert rep["injection_order"]["status"] == "biased" and "may_be_file_order" not in rep["injection_order"]


def test_compare_confounding_reads_the_filename_declared_source():
    """比較の交絡の判定（compare が組む代理の単位）にも、単位の注入順の出所が届く。"""
    from dataclasses import replace

    from metabolomix.analysis.stats.compare import compare
    from tests.test_stats_preprocess import _drift_study

    study = _drift_study()
    study.samples = [replace(s, injection_order_source="filename") for s in study.samples]
    cfg = {"filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None},
           "normalize": {"method": "none"}, "drift": {"mode": "off"}}
    prep = P.preprocess(study, cfg)
    groups = sorted({u.group for u in prep.units if u.group is not None})
    res = compare(prep, factors=["group"],
                  comparisons=[{"id": "x", "factor": "group", "reference": groups[0], "test": groups[1]}])
    assert res["confounding"]["injection_order"]["may_be_file_order"] is True
    assert any(w["code"] == "injection_order_may_be_file_order" for w in res["warnings"])
