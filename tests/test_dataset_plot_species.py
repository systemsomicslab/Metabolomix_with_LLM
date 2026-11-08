"""dataset_plot の group_intensity・species・species_pca（B-2b Task 7）。

旧 `test_group_intensity_tools.py` / `test_species_tools.py`（ARF の行から選ぶ旧ツール。Task 8 で消える）の
うち、データの意味に関わる主張——キュレーションの当て方・検出限界・MS/MS の要求・share の分母・標準液・
除外・gap-fill の印——を `DatasetState`（解析フォルダを `load_dataset` で開いたもの）の上へ移したもの。
fixture はテスト自身が tmp に作る（`.arf` `.arf2` と param）。
"""
from __future__ import annotations

import json

import numpy as np
import pytest
from mcp.server.fastmcp import Image

from metabolomix.core import handles, session_state
from metabolomix.curation import flags as curation_flags
from metabolomix.dataset import msdial_adapter
from metabolomix.plots import dataset_rows
from metabolomix.tools import handle_tools, plot_tools, reports
from tests.config_fixtures import use_config
from tests.curation_fixtures import arf2_spot_raw, arf_row, match_result, write_arf, write_arf2

STEM = "AlignmentResult_2026_01_01_00_00_00"
SAMPLES = ["blank_1", "ctrl_1", "ctrl_2", "ctrl_3", "ko_1", "ko_2", "ko_3", "std_1"]
NAME = {n: f"20261006_{n}" for n in SAMPLES}
PNG_MAGIC = b"\x89PNG"
MSMS = {1: 1}      # `.arf` の行の Key10（MS2RawSpectrumID2CE）。非空なら MS/MS を取得したセル


def _group(heights, gap=(), msms=False):
    rows = []
    for fid, (name, h) in enumerate(zip(SAMPLES, heights)):
        row = arf_row(file_id=fid, mz=700.0, rt=3.0, height=h, gap_filled=name in gap)
        row[1] = NAME[name]
        if msms:
            row[10] = dict(MSMS)
        rows.append(row)
    return rows


#            blank ctrl_1 ctrl_2 ctrl_3 ko_1 ko_2 ko_3 std_1
SPOTS = [
    # (spot_id 0..) 高さ、arf2 の Name / Ontology / adduct、MS/MS の有無、gap-fill の試料
    ([10, 5000, 6000, 5500, 2000, 2500, 2200, 10], "PG 34:1|PG 16:0_18:1", "PG", "[M-H]-", True, ("ctrl_2",)),
    ([10, 3000, 3500, 3200, 4000, 4200, 4100, 10], "PG 35:1|PG 16:0_19:1", "PG", "[M-H]-", False, ()),
    ([1500, 1400, 1300, 1450, 1500, 1500, 1500, 10], "PC 33:1(d7)|PC 15:0_18:1(d7)", "PC", "[M+H]+", True, ()),
    ([5, 800, 900, 700, 600, 650, 620, 5], "no MS2: PC 34:1", "PC", "[M+H]+", True, ()),
    ([0, 20, 30, 25, 20, 22, 21, 9000], "PG 31:1|PG 17:0_14:1", "PG", "[M-H]-", False, ()),   # 標準液だけ
    ([0, 2000, 2100, 1900, 4000, 4300, 300, 0], "DGDG 34:1|DGDG 16:0_18:1", "DGDG", "[M+CH3COO]-", True, ()),
    ([0, 900, 950, 920, 800, 700, 760, 0], "DG 34:1|DG 16:0_18:1", "DG", "[M+NH4]+", True, ()),
    ([0, 300, 320, 310, 280, 260, 270, 0], "DG 34:1|DG 16:0_18:1", "DG", "[M+Na]+", True, ()),
]


@pytest.fixture(autouse=True)
def fresh_state(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": tmp_path / "reports"})
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


@pytest.fixture()
def loaded(tmp_path):
    folder = tmp_path / "msdial"
    folder.mkdir()
    write_arf(folder / f"{STEM}_PeakProperties.arf",
              [_group(h, gap=gap, msms=msms) for h, _n, _o, _a, msms, gap in SPOTS])
    arf2 = write_arf2(folder / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=i, name=name, ontology=onto, adduct=adduct, mz=700.0, rt=3.0,
                      matches=[match_result()] if msms else [])
        for i, (_h, name, onto, adduct, msms, _g) in enumerate(SPOTS)])
    body = json.loads(handle_tools.load_dataset(directory=str(folder)))
    assert body["status"] == "success", body
    handle = handles.registry.get(body["dataset_handle"])
    return {"ds": handle.dataset, "arf2": arf2, "folder": folder, "handle": handle}


def _p(kind, **kw):
    out = plot_tools.dataset_plot(kind=kind, output="payload", **kw)
    assert isinstance(out, str), out
    return json.loads(out)


def _gi(**kw):
    kw.setdefault("groups", ["ctrl", "ko"])
    return _p("group_intensity", **kw)


def _species(**kw):
    kw.setdefault("groups", ["ctrl", "ko"])
    return _p("species", **kw)


def _pca(**kw):
    kw.setdefault("groups", ["ctrl", "ko"])
    return _p("species_pca", **kw)


def _ids(item):
    return [s["spot_id"] for s in item["spots"]]


def _error(out):
    return json.loads(out)["error"]


def _write_flags(loaded, entries):
    store = curation_flags.FlagStore(curation_flags.curation_dir(loaded["arf2"]))
    store.append(entries, alignment=curation_flags.alignment_key(loaded["arf2"]),
                 review_id=None, source="test")
    msdial_adapter.reload_curation(loaded["ds"])


def _write_review(arf2, review_id, sha, spots):
    directory = curation_flags.curation_dir(arf2)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"review-{review_id}.json").write_text(
        json.dumps({"alignment": {"alignment_sha256": sha}, "spots": spots}), encoding="utf-8")


# ---------- group_intensity: 項目の解決・除外 ----------

def test_missing_dataset_returns_missing_state():
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"]))
    assert err["code"] == "missing_state" and err["required_tools"] == ["load_dataset"]


def test_payload_classes_names_ids_and_nd(loaded):
    p = _gi(items=["PG", "PC", "PE", "PG 34:1", "1"])
    assert p["plot_schema"] == "lipidmix.group_intensity.v1"
    by = {it["item"]: it for it in p["items"]}
    assert _ids(by["PG"]) == [0, 1, 4]
    assert _ids(by["PC"]) == [3]                       # 内部標準は除く
    assert by["PC"]["n_spots_msms"] == 0               # `no MS2:` は MS/MS の裏付けに数えない
    assert by["PG"]["n_spots_msms"] == 1
    assert by["PE"]["detected"] is False
    assert _ids(by["PG 34:1"]) == [0]
    assert _ids(by["1"]) == [1] and by["1"]["parts"][0]["kind"] == "feature"      # 特徴 ID はそのまま
    assert [g["label"] for g in p["groups"]] == ["ctrl", "ko"]
    assert p["groups"][0]["samples"] == [NAME["ctrl_1"], NAME["ctrl_2"], NAME["ctrl_3"]]
    assert p["detection_limit"] is None


def test_values_are_the_class_sum_of_the_matrix(loaded):
    p = _gi(items=["PG"], standard_samples=["std"])
    ctrl = {s["sample"]: s["value"] for s in p["items"][0]["groups"][0]["samples"]}
    assert ctrl[NAME["ctrl_1"]] == 5000 + 3000 and ctrl[NAME["ctrl_3"]] == 5500 + 3200


def test_standard_samples_drop_standard_only_spots(loaded):
    p = _gi(items=["PG"], standard_samples=["std"])
    assert _ids(p["items"][0]) == [0, 1]
    assert p["excluded"]["standard_only"] == ["#4 PG 31:1|PG 17:0_14:1"]


def test_standard_filter_still_works_when_standard_sample_is_excluded(loaded):
    loaded["ds"].excluded_samples = {NAME["std_1"]: "x"}
    p = _gi(items=["PG"], standard_samples=["std"])
    assert _ids(p["items"][0]) == [0, 1]
    assert p["excluded"]["standard_only"] == ["#4 PG 31:1|PG 17:0_14:1"]


def test_curation_flags_are_applied_by_default(loaded):
    _write_flags(loaded, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    p = _gi(items=["PG"], groups=["ctrl"])
    assert 1 not in _ids(p["items"][0]) and p["excluded"]["curation"] == ["#1 PG 35:1|PG 16:0_19:1"]
    p = _gi(items=["PG"], groups=["ctrl"], apply_curation=False)
    assert 1 in _ids(p["items"][0])


def test_curation_assign_moves_a_spot_to_another_class(loaded):
    _write_flags(loaded, [{"spot_id": 1, "flag": "assign", "name": "PE 34:1", "ontology": "PE",
                           "inchikey": "", "note": ""}])
    p = _gi(items=["PG", "PE"], groups=["ctrl"])
    by = {it["item"]: it for it in p["items"]}
    assert 1 not in _ids(by["PG"]) and _ids(by["PE"]) == [1]


def test_auto_likely_wrong_with_review_excludes_spots_and_ignores_other_alignment(loaded):
    arf2 = loaded["arf2"]
    sha = curation_flags.alignment_key(arf2)["alignment_sha256"]
    _write_review(arf2, "cr-20261006-120000-aaaa", sha, [{"spot_id": 1, "verdict": "likely_wrong"}])
    # 別アラインメント（sha256 違い）のレビューは、より新しくても無視される
    _write_review(arf2, "cr-20261006-130000-bbbb", "0" * 64, [{"spot_id": 0, "verdict": "likely_wrong"}])
    p = _gi(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert _ids(p["items"][0]) == [0, 4]
    assert p["excluded"]["auto_likely_wrong"][0].startswith("#1 ")
    assert not any("レビュー" in c for c in p["caveats"])
    assert 1 in _ids(_gi(items=["PG"], groups=["ctrl"])["items"][0])


def test_auto_likely_wrong_without_review_adds_caveat(loaded):
    p = _gi(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert any("レビュー" in c for c in p["caveats"])


def test_auto_likely_wrong_keeps_confirmed_spots(loaded):
    from metabolomix.msdial import tags as msdial_tags

    arf2 = loaded["arf2"]
    sha = curation_flags.alignment_key(arf2)["alignment_sha256"]
    _write_review(arf2, "cr-20261006-120000-aaaa", sha, [{"spot_id": 1, "verdict": "likely_wrong"}])
    msdial_tags.update_alignment_tag(msdial_tags.alignment_tag_path(arf2),
                                     tag_id=msdial_tags.CONFIRMED_TAG_ID, add=[1], remove=[])
    p = _gi(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert 1 in _ids(p["items"][0]) and p["excluded"]["auto_likely_wrong"] == []


def test_broken_flag_file_returns_error_with_message(loaded):
    """人の判断は読み込みのとき（`ds.curation_flags`）に読む。壊れた記録は、最新レビューを読み直す
    `exclude_auto_likely_wrong` のところで名指しで止まる。"""
    arf2 = loaded["arf2"]
    sha = curation_flags.alignment_key(arf2)["alignment_sha256"]
    _write_review(arf2, "cr-20261006-120000-aaaa", sha, [{"spot_id": 1, "verdict": "likely_wrong"}])
    (curation_flags.curation_dir(arf2) / curation_flags.FLAGS_FILENAME).write_text("{not json\n", encoding="utf-8")
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"],
                                         exclude_auto_likely_wrong=True))
    assert err["code"] == "DATASET_BAD_REQUEST" and err["message"] and err["details"]["flags_file"]


def test_broken_flag_file_at_reload_keeps_the_previous_judgements(loaded):
    _write_flags(loaded, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    (curation_flags.curation_dir(loaded["arf2"]) / curation_flags.FLAGS_FILENAME).write_text(
        "{not json\n", encoding="utf-8")
    msdial_adapter.reload_curation(loaded["ds"])
    assert any("キュレーションのフラグを読めなかった" in n for n in loaded["ds"].load_notes)
    assert 1 not in _ids(_gi(items=["PG"], groups=["ctrl"])["items"][0])      # 前の判断が残る


def test_apply_curation_without_a_record_adds_caveat(loaded):
    loaded["ds"].curation_flags = None
    p = _gi(items=["PG"], groups=["ctrl"])
    assert any("キュレーションの記録" in c for c in p["caveats"])


def test_excluded_samples_and_low_reliability(loaded):
    loaded["ds"].excluded_samples = {NAME["ko_2"]: "x"}
    p = _gi(items=["PG"], low_reliability_samples=["ctrl_2"])
    assert p["groups"][1]["samples"] == [NAME["ko_1"], NAME["ko_3"]]
    ctrl = p["items"][0]["groups"][0]
    assert ctrl["n_in_stats"] == 2
    assert [s["low_reliability"] for s in ctrl["samples"]] == [False, True, False]


def test_group_fully_excluded_is_an_error(loaded):
    loaded["ds"].excluded_samples = {NAME[n]: "x" for n in ("ko_1", "ko_2", "ko_3")}
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl", "ko"]))
    assert err["code"] == "DATASET_BAD_REQUEST" and "No sample matched spec 'ko'" in err["message"]


def test_sample_metadata_exclusion_include_false_is_not_drawn(loaded):
    ds = loaded["ds"]
    rows = ds.sample_metadata_rows or __import__(
        "metabolomix.analysis.sample_manifest", fromlist=["x"]).resolve_metadata(ds, None)
    ds.sample_metadata_rows = rows
    next(r for r in rows if r["sample_id"] == NAME["ko_3"])["include"] = False
    p = _gi(items=["PG"])
    assert NAME["ko_3"] not in p["groups"][1]["samples"]


def test_excluded_features_are_not_drawn_in_class_or_name_or_id_parts(loaded):
    loaded["ds"].excluded_features = {"1": "x"}
    p = _gi(items=["PG", "PG 35:1", "1"], groups=["ctrl"])
    by = {it["item"]: it for it in p["items"]}
    assert 1 not in _ids(by["PG"])
    assert by["PG 35:1"]["spots"] == [] and by["PG 35:1"]["detected"] is False
    assert by["1"]["spots"] == []
    assert p["excluded"]["manual"] == ["#1 PG 35:1|PG 16:0_19:1"]


def test_caption_reports_manual_exclusions(loaded):
    loaded["ds"].excluded_features = {"1": "x"}
    out = plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"])
    assert "手動除外 1" in out[0]


def test_detection_limit_from_param_file_and_argument(loaded):
    (loaded["folder"] / "Dataset_x_param_202610061200.txt").write_text(
        "Minimum peak height: 1000\n", encoding="utf-8")
    p = _gi(items=["PG"], groups=["ctrl"])
    assert p["detection_limit"] == 1000.0 and p["detection_limit_source"] == "param_file"
    p = _gi(items=["PG"], groups=["ctrl"], detection_limit=500)
    assert p["detection_limit"] == 500.0 and p["detection_limit_source"] == "argument"


@pytest.mark.parametrize("limit", [-5, 0, float("nan"), float("inf")])
def test_invalid_detection_limit_is_an_error(loaded, limit):
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"],
                                         detection_limit=limit))
    assert err["code"] == "DATASET_BAD_REQUEST" and "detection_limit" in err["message"]


def test_gap_filled_cells_marked_in_group_intensity(loaded):
    """ctrl_2 の PG 34:1（spot 0）は gap-filled。検出マスクの False が印になる。"""
    ds = loaded["ds"]
    assert ds.detected_mask is not None and not ds.detected_mask[0][SAMPLES.index("ctrl_2")]
    p = _gi(items=["PG 34:1", "PG 35:1"], groups=["ctrl"])
    by = {it["item"]: it for it in p["items"]}
    gap = {s["sample"]: s["gap_filled_fraction"] for s in by["PG 34:1"]["groups"][0]["samples"]}
    assert gap[NAME["ctrl_2"]] == 1.0 and gap[NAME["ctrl_1"]] == 0.0
    assert all(s["gap_filled_fraction"] == 0.0 for s in by["PG 35:1"]["groups"][0]["samples"])


def test_missing_detected_mask_adds_caveat_and_marks_nothing(loaded):
    loaded["ds"].detected_mask = None
    p = _gi(items=["PG 34:1"], groups=["ctrl"])
    assert any("detected_mask" in c for c in p["caveats"])
    assert all(s["gap_filled_fraction"] == 0.0 for s in p["items"][0]["groups"][0]["samples"])


@pytest.mark.parametrize("kw", [{"items": []}, {"groups": []}, {"groups": ["nonexistent"]},
                                {"items": [" + "]}])
def test_invalid_input_is_bad_request(loaded, kw):
    args = {"items": ["PG"], "groups": ["ctrl"], **kw}
    err = _error(plot_tools.dataset_plot(kind="group_intensity", **args))
    assert err["code"] == "DATASET_BAD_REQUEST" and err["message"]


# ---------- species ----------

def test_share_payload_expands_classes_and_uses_items_as_default_basis(loaded):
    p = _species(items=["PG", "DGDG"], standard_samples=["std"])
    assert p["plot_schema"] == "lipidmix.species_intensity.v1" and p["value"] == "share"
    assert [s["spot_id"] for s in p["spots"]] == [0, 1, 5]
    assert p["share_basis"]["spot_ids"] == [0, 1, 5]
    c1 = p["spots"][0]["groups"][0]["samples"][0]
    assert c1["share"] == pytest.approx(5000 / (5000 + 3000 + 2000) * 100, abs=1e-3)


def test_share_denominator_is_the_basis_total_and_changes_with_share_basis(loaded):
    p = _species(items=["PG"], standard_samples=["std"], share_basis=["PG", "DGDG"])
    assert p["share_basis"]["spot_ids"] == [0, 1, 5]
    c1 = p["spots"][0]["groups"][0]["samples"][0]
    assert c1["share"] == pytest.approx(5000 / 10000 * 100, abs=1e-3)
    assert [s["spot_id"] for s in p["spots"]] == [0, 1]


def test_share_basis_and_require_msms(loaded):
    p = _species(items=["PG"], standard_samples=["std"], share_basis=["PG", "DGDG"], require_msms=True)
    assert [s["spot_id"] for s in p["spots"]] == [0]
    assert p["share_basis"]["spot_ids"] == [0, 5]        # 分母にも MS/MS の条件が効く
    assert p["excluded"]["no_msms"] == ["#1 PG 35:1|PG 16:0_19:1"]
    assert p["spots"][0]["msms"] is True


def test_require_msms_without_the_evidence_column_is_explained(loaded):
    ds = loaded["ds"]
    for meta in ds.feature_metadata.values():
        meta.pop("msms_assigned", None)
    err = _error(plot_tools.dataset_plot(kind="species", items=["PG"], groups=["ctrl"], require_msms=True))
    assert err["code"] == "DATASET_BAD_REQUEST" and "require_msms" in err["message"]
    assert _species(items=["PG"], groups=["ctrl"])["spots"]        # 要求しなければ描ける


def test_basis_without_the_plotted_species_adds_caveat(loaded):
    p = _species(items=["PG"], share_basis=["DGDG"])
    assert any("分母に含まれない" in c for c in p["caveats"])


def test_unmatched_share_basis_part_is_an_error(loaded):
    err = _error(plot_tools.dataset_plot(kind="species", items=["PG", "PE"], groups=["ctrl", "ko"],
                                         share_basis=["PG", "DGDGX"]))
    assert err["code"] == "DATASET_BAD_REQUEST" and "DGDGX" in err["message"] and "share_basis" in err["message"]
    assert "species" not in loaded["ds"].figures


def test_unmatched_item_part_is_a_caveat_in_payload_and_caption(loaded):
    p = _species(items=["PG", "NOSUCH"])
    assert any("NOSUCH" in c for c in p["caveats"])
    out = plot_tools.dataset_plot(kind="species", items=["PG", "NOSUCH"], groups=["ctrl", "ko"])
    assert "NOSUCH" in out[0]


def test_intensity_with_share_basis_warns_that_basis_is_ignored(loaded):
    p = _species(items=["PG"], value="intensity", share_basis=["PG", "DGDG"])
    assert p["value"] == "height"                    # 形式（lipidmix.species_intensity.v1）の値の名前
    assert any("share_basis" in c and "無視" in c for c in p["caveats"])


def test_same_name_different_adduct_are_separate_panels(loaded):
    p = _species(items=["DG"], value="intensity")
    assert [(s["spot_id"], s["adduct"]) for s in p["spots"]] == [(6, "[M+NH4]+"), (7, "[M+Na]+")]


def test_species_gap_filled_cells_are_marked_and_low_reliability_excluded_from_stats(loaded):
    p = _species(items=["PG 34:1"], low_reliability_samples=["ko_3"])
    ctrl = p["spots"][0]["groups"][0]["samples"]
    assert [s["gap_filled"] for s in ctrl] == [False, True, False]
    ko = p["spots"][0]["groups"][1]
    assert ko["n_in_stats"] == 2 and ko["samples"][2]["low_reliability"] is True


def test_too_many_panels_is_an_error_and_keeps_the_previous_figure(loaded, monkeypatch):
    from metabolomix.plots import species

    monkeypatch.setattr(species, "MAX_PANELS", 2)
    plot_tools.dataset_plot(kind="species", items=["PG 34:1"], groups=["ctrl"])
    before = loaded["ds"].figures["species"]
    err = _error(plot_tools.dataset_plot(kind="species", items=["PG", "DGDG"], groups=["ctrl"],
                                         standard_samples=["std"]))
    assert err["code"] == "DATASET_BAD_REQUEST" and "40" not in err["message"] and "2" in err["message"]
    assert loaded["ds"].figures["species"] is before


@pytest.mark.parametrize("kw", [{"value": "area"}, {"value": "height"}, {"groups": ["nonexistent"]},
                                {"items": ["NOSUCH"]}])
def test_species_invalid_input_is_bad_request(loaded, kw):
    err = _error(plot_tools.dataset_plot(kind="species", **{"items": ["PG"], "groups": ["ctrl"], **kw}))
    assert err["code"] == "DATASET_BAD_REQUEST" and err["message"]


def test_species_image_caption_and_figure_record(loaded):
    out = plot_tools.dataset_plot(kind="species", items=["PG"], groups=["ctrl", "ko"],
                                  low_reliability_samples=["ko_3"], standard_samples=["std"])
    assert isinstance(out, list) and isinstance(out[1], Image) and out[1].data[:4] == PNG_MAGIC
    assert "2 パネル" in out[0] and "ctrl n=3" in out[0]
    record = loaded["ds"].figures["species"]
    assert record["payload"]["spots"][0]["spot_id"] == 0
    assert record["metadata_revision"] == loaded["ds"].metadata_revision


# ---------- species_pca ----------

def test_species_pca_projects_low_reliability_and_orients(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], low_reliability_samples=["ko_3"], orient_by="ko",
             standard_samples=["std"])
    assert p["plot_schema"] == "lipidmix.species_pca.v1"
    pts = {pt["label"]: pt for pt in p["points"]}
    assert pts[NAME["ko_3"]]["fitted"] is False and pts[NAME["ctrl_1"]]["fitted"] is True
    ko_x = [pts[NAME[f"ko_{i}"]]["x"] for i in (1, 2)]
    assert sum(ko_x) / 2 > 0                                    # orient_by の群が正
    assert p["n_fit"] == 5 and p["scaling"] == "autoscale" and p["n_species"] == 5
    assert "loadings_rows" not in p                              # 全行は結果に持つ（payload を肥大させない）
    result = loaded["ds"].results[p["result_id"]]
    assert {row["feature_id"] for row in result["loadings_rows"]} == {"0", "1", "5", "6", "7"}
    assert all(-1.0001 <= v <= 1.0001 for row in result["loadings_rows"] for v in row["r"])
    assert result["provenance"]["kind"] == "species_pca"


def test_species_pca_deterministic_over_500_features(loaded, monkeypatch):
    """分子種 PCA は 500 超の特徴でも randomized SVD に落ちないよう `svd_solver="full"` を渡す。"""
    from metabolomix.tools import plot_species_kinds

    seen = []
    real = plot_species_kinds.run_pca_fit_subset

    def spy(*args, **kwargs):
        seen.append(kwargs.get("svd_solver"))
        return real(*args, **kwargs)

    monkeypatch.setattr(plot_species_kinds, "run_pca_fit_subset", spy)
    _pca(items=["PG", "DGDG", "DG"])
    assert seen == ["full"]


def test_species_pca_registers_a_result_and_a_figure(loaded):
    ds = loaded["ds"]
    p = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])
    record = ds.figures["species_pca"]
    assert record["result_id"] == p["result_id"] and record["result"] is ds.results[p["result_id"]]
    assert record["metadata_revision"] == ds.metadata_revision
    # スコア図の kind（pca）には使えない。loadings は受ける
    err = _error(plot_tools.dataset_plot(kind="pca", result_id=p["result_id"]))
    assert err["code"] == "DATASET_BAD_REQUEST"


def test_species_pca_loadings_are_drawable_and_savable(loaded, tmp_path):
    p = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])
    rid = p["result_id"]
    out = _p("loadings", result_id=rid, top_n=3)
    assert out["plot_schema"] == "lipidmix.pca_loadings.v1" and out["result_id"] == rid
    labels = {r["label"] for panel in out["panels"] for r in panel["rows"]}
    assert labels <= {"PG 16:0_18:1", "PG 16:0_19:1", "DGDG 16:0_18:1", "DG 16:0_18:1"}
    image = plot_tools.dataset_plot(kind="loadings", result_id=rid)
    assert isinstance(image[1], Image)
    assert "l-1_loadings.png" in reports.save_figure("loadings", "l-1", result_id=rid)
    assert (tmp_path / "reports" / "figures" / "l-1_loadings.svg").is_file()
    # result_id を省いた loadings は直近の dataset_pca を使う（分子種 PCA の結果は使わない）
    err = _error(plot_tools.dataset_plot(kind="loadings"))
    assert err["code"] == "missing_state"


def test_species_pca_result_goes_stale_after_exclusion_or_metadata_change(loaded, tmp_path):
    from metabolomix.tools.sample_tools import dataset_exclude

    ds = loaded["ds"]
    rid = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])["result_id"]
    plot_tools.dataset_plot(kind="loadings", result_id=rid)
    body = json.loads(dataset_exclude(samples=[NAME["ko_3"]], reason="test"))
    assert body.get("changed") is True, body
    for call in (lambda: plot_tools.dataset_plot(kind="loadings", result_id=rid),
                 lambda: reports.save_figure("loadings", "s-1", result_id=rid)):
        assert _error(call())["code"] == "STALE_ANALYSIS_RESULT"
    # キュレーションの判断の更新では古くしない
    rid2 = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])["result_id"]
    _write_flags(loaded, [{"spot_id": 1, "flag": "wrong", "note": ""}])
    assert _p("loadings", result_id=rid2)["result_id"] == rid2
    ds.metadata_revision += 1
    assert _error(plot_tools.dataset_plot(kind="loadings", result_id=rid2))["code"] == "STALE_ANALYSIS_RESULT"


def test_species_pca_total_normalization_with_zero_total_names_the_sample(loaded):
    err = _error(plot_tools.dataset_plot(kind="species_pca", items=["DG"], groups=["blank", "ctrl"],
                                         normalize="total"))
    assert err["code"] == "DATASET_BAD_REQUEST" and "blank_1" in err["message"]


def test_species_pca_sample_in_two_groups_is_counted_once(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], groups=["ctrl", "ctrl_1", "ko"], standard_samples=["std"])
    labels = [pt["label"] for pt in p["points"]]
    assert labels.count(NAME["ctrl_1"]) == 1


@pytest.mark.parametrize("kw", [{"normalize": "median"}, {"orient_by": "nope"}, {"items": ["PG 34:1"]},
                                {"groups": ["ctrl"], "low_reliability_samples": ["ctrl_1", "ctrl_2"]},
                                {"orient_by": "ko", "low_reliability_samples": ["ko"]}])
def test_species_pca_invalid_input_is_bad_request_and_keeps_the_previous_figure(loaded, kw):
    _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])
    before = loaded["ds"].figures["species_pca"]
    n_results = len(loaded["ds"].results)
    args = {"items": ["PG", "DGDG", "DG"], "groups": ["ctrl", "ko"], "standard_samples": ["std"], **kw}
    err = _error(plot_tools.dataset_plot(kind="species_pca", **args))
    assert err["code"] == "DATASET_BAD_REQUEST" and err["message"]
    assert loaded["ds"].figures["species_pca"] is before and len(loaded["ds"].results) == n_results


def test_species_pca_image_and_caption(loaded):
    out = plot_tools.dataset_plot(kind="species_pca", items=["PG", "DGDG", "DG"], groups=["ctrl", "ko"],
                                  standard_samples=["std"])
    assert isinstance(out[1], Image) and out[1].data[:4] == PNG_MAGIC and "PC1" in out[0]


def test_species_pca_caption_carries_review_caveat_and_exclusion_summary(loaded):
    out = plot_tools.dataset_plot(kind="species_pca", items=["PG", "DGDG", "DG"], groups=["ctrl", "ko"],
                                  standard_samples=["std"], exclude_auto_likely_wrong=True)
    assert "curation_review" in out[0] and "除外 — " in out[0] and "注意:" in out[0]


def test_species_pca_settings_record_filters(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], exclude_auto_likely_wrong=True, apply_curation=False,
             standard_samples=["std"])
    result = loaded["ds"].results[p["result_id"]]
    for settings in (p["settings"], result["provenance"]["effective_parameters"]):
        assert settings["apply_curation"] is False and settings["exclude_auto_likely_wrong"] is True
        assert settings["standard_samples"] == ["std"]


def test_species_pca_top_and_bottom_do_not_overlap(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])
    top, bottom = p["pc1_r_top"], p["pc1_r_bottom"]
    assert not set(top) & set(bottom)
    assert all(" +" in t or "+0" in t for t in top) and all(" -" in b for b in bottom)
    assert len(top) + len(bottom) <= 5 * 2


def test_species_pca_respects_excluded_samples_and_features(loaded):
    ds = loaded["ds"]
    ds.excluded_samples = {NAME["ko_3"]: "x"}
    ds.excluded_features = {"5": "x"}
    p = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])
    assert NAME["ko_3"] not in {pt["label"] for pt in p["points"]}
    assert p["n_species"] == 4 and p["excluded"]["manual"] == ["#5 DGDG 34:1|DGDG 16:0_18:1"]


# ---------- 引数の持ち主・形式名・図の古さ・save_figure ----------

def test_argument_ownership_names_the_owning_kinds(loaded):
    err = _error(plot_tools.dataset_plot(kind="species", items=["PG"], groups=["ctrl"], detection_limit=5))
    assert err["code"] == "DATASET_BAD_REQUEST" and "group_intensity" in err["message"]
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"],
                                         share_basis=["PG"]))
    assert "species" in err["message"] and err["details"]["belongs_to"] == ["species"]
    err = _error(plot_tools.dataset_plot(kind="species", items=["PG"], groups=["ctrl"], normalize="total"))
    assert "species_pca" in err["message"]
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"],
                                         require_msms=True))
    assert err["details"]["belongs_to"] == ["species", "species_pca"]       # 共有する kind は全部挙げる
    err = _error(plot_tools.dataset_plot(kind="volcano", items=["PG"]))
    assert err["details"]["belongs_to"] == ["group_intensity", "species", "species_pca"]
    # value は loadings と species が共有する（loadings の値は r / coefficient、species の値は share / intensity）
    err = _error(plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"], value="share"))
    assert err["details"]["belongs_to"] == ["loadings", "species"]


def test_figures_with_baked_groups_go_stale_when_the_metadata_revision_moves(loaded, tmp_path):
    ds = loaded["ds"]
    for kind, args in (("group_intensity", {}), ("species", {}),
                       ("species_pca", {"items": ["PG", "DGDG", "DG"], "standard_samples": ["std"]})):
        plot_tools.dataset_plot(kind=kind, **{"items": ["PG"], "groups": ["ctrl", "ko"], **args})
        assert "a-1_" + kind + ".png" in reports.save_figure(kind, "a-1")
    figures = tmp_path / "reports" / "figures"
    assert (figures / "a-1_group_intensity.svg").is_file() and (figures / "a-1_species_pca.svg").is_file()
    ds.metadata_revision += 1                               # 試料の情報の更新
    for kind in ("group_intensity", "species", "species_pca"):
        err = _error(reports.save_figure(kind, "b-1"))
        assert err["code"] == "missing_state" and "dataset_plot" in err["required_tools"]
        assert "作り直" in err["message"]
        assert not (figures / f"b-1_{kind}.png").exists()


def test_species_pca_figure_is_dropped_with_the_pca_set_and_save_figure_rechecks(loaded):
    from metabolomix.analysis import result_state

    ds = loaded["ds"]
    p = _pca(items=["PG", "DGDG", "DG"], standard_samples=["std"])
    assert "species_pca" in ds.figures
    # 取りこぼした古さ（結果が現在でない）は save_figure の直前の確認が塞ぐ
    ds.results[p["result_id"]]["provenance"]["dataset_id"] = "ds_other"
    err = _error(reports.save_figure("species_pca", "c-1"))
    assert err["code"] == "missing_state" and "作り直" in err["message"]
    ds.results[p["result_id"]]["provenance"]["dataset_id"] = ds.dataset_id
    # 除外（前処理の入力の変更）で図ごと捨てられる
    result_state.invalidate_results(ds, {"exclusion"})
    assert "species_pca" not in ds.figures
    err = _error(reports.save_figure("species_pca", "c-1"))
    assert err["code"] == "missing_state"


def test_save_figure_result_id_is_not_accepted_for_the_new_kinds(loaded):
    err = json.loads(reports.save_figure("species", "x", result_id="r1"))
    assert err["status"] == "error" and "result_id" in err["message"]


def test_failed_render_does_not_store_the_figure(loaded, monkeypatch):
    from metabolomix.plots import group_intensity as gi

    def boom(*a, **k):
        raise RuntimeError("render failed")
    monkeypatch.setattr(gi, "render_group_intensity_plot", boom)
    with pytest.raises(RuntimeError):
        plot_tools.dataset_plot(kind="group_intensity", items=["PG"], groups=["ctrl"])
    assert "group_intensity" not in loaded["ds"].figures


# ---------- dataset_rows（部品） ----------

def test_select_items_expands_classes_and_passes_ids_through(loaded):
    ds = loaded["ds"]
    assert dataset_rows.select_items(ds, ["DG"], apply_curation=True, exclude_auto_likely_wrong=False,
                                     require_msms=False) == ["6", "7"]
    assert dataset_rows.select_items(ds, ["5", "PG 34:1"], apply_curation=True,
                                     exclude_auto_likely_wrong=False, require_msms=False) == ["5", "0"] or \
        sorted(dataset_rows.select_items(ds, ["5", "PG 34:1"], apply_curation=True,
                                         exclude_auto_likely_wrong=False, require_msms=False)) == ["0", "5"]
    ds.excluded_features = {"6": "x"}
    assert dataset_rows.select_items(ds, ["DG"], apply_curation=True, exclude_auto_likely_wrong=False,
                                     require_msms=False) == ["7"]
    assert dataset_rows.select_items(ds, ["PG"], apply_curation=True, exclude_auto_likely_wrong=False,
                                     require_msms=True) == ["0"]


def test_rows_by_spot_shape_values_and_exclusions(loaded):
    ds = loaded["ds"]
    ds.excluded_samples = {NAME["ctrl_3"]: "x"}
    rows = dataset_rows.rows_by_spot(ds, ["0"], groups=["ctrl", "ko"])
    by = {r["file_name"]: r for r in rows[0]}
    assert NAME["ctrl_3"] not in by and NAME["std_1"] not in by          # 除外した試料・描かない試料は行にしない
    assert by[NAME["ctrl_1"]] == {"file_name": NAME["ctrl_1"], "height": 5000.0, "is_gap_filled": False}
    assert by[NAME["ctrl_2"]]["is_gap_filled"] is True


def test_rows_by_spot_skips_nan_cells(loaded):
    ds = loaded["ds"]
    ds.feature_matrix = np.array(ds.feature_matrix, dtype=float)
    ds.feature_matrix[0][SAMPLES.index("ctrl_1")] = np.nan
    rows = dataset_rows.rows_by_spot(ds, ["0"], groups=["ctrl"])
    assert NAME["ctrl_1"] not in {r["file_name"] for r in rows[0]}
