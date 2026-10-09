# tests/test_group_intensity_tools.py
"""arf_plot_group_intensity / save_figure(kind="group_intensity")（ツール層）。fixture はテストが作る。"""
from __future__ import annotations

import json

import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state
from metabolomix.curation import flags as curation_flags
from tests.curation_fixtures import arf2_spot_raw, arf_row, match_result, write_arf, write_arf2

STEM = "AlignmentResult_2026_01_01_00_00_00"
SAMPLES = ["blank_1", "ctrl_1", "ctrl_2", "ko_1", "ko_2", "std_1"]


def _group(heights, gap=()):
    rows = []
    for fid, (name, h) in enumerate(zip(SAMPLES, heights)):
        row = arf_row(file_id=fid, mz=700.0, rt=3.0, height=h, gap_filled=name in gap)
        row[1] = f"20261006_{name}"
        rows.append(row)
    return rows


@pytest.fixture()
def loaded(tmp_path):
    arf = write_arf(tmp_path / f"{STEM}_PeakProperties.arf", [
        _group([10, 5000, 6000, 2000, 2500, 10]),            # 0 PG 34:1
        _group([10, 3000, 3500, 4000, 4200, 10]),            # 1 PG 35:1
        _group([1500, 1400, 1300, 1450, 1500, 10]),          # 2 PC 33:1(d7) 内部標準
        _group([5, 800, 900, 700, 600, 5]),                  # 3 no MS2: PC 34:1
        _group([0, 20, 30, 25, 20, 9000]),                   # 4 PG 31:1 標準液だけ
    ])
    arf2 = write_arf2(tmp_path / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=0, name="PG 34:1|PG 16:0_18:1", ontology="PG"),
        arf2_spot_raw(spot_id=1, name="PG 35:1|PG 16:0_19:1", ontology="PG"),
        arf2_spot_raw(spot_id=2, name="PC 33:1(d7)|PC 15:0_18:1(d7)", ontology="PC"),
        arf2_spot_raw(spot_id=3, name="no MS2: PC 34:1", ontology="PC",
                      matches=[match_result({3: -1.0})]),
        arf2_spot_raw(spot_id=4, name="PG 31:1|PG 17:0_14:1", ontology="PG"),
    ])
    session_state.session = session_state.AnalysisSession()
    with open(arf, "rb") as handle:
        session_state.session.arf.features = arf_reader.deserialize(handle)
    session_state.session.arf.current_file_path = str(arf)
    session_state.session.arf.class_index = None
    return {"arf": arf, "arf2": arf2}


def _payload(**kw):
    from metabolomix.arf.tools import arf_plot_group_intensity
    return json.loads(arf_plot_group_intensity(output="payload", **kw))


def test_missing_arf_returns_missing_state():
    from metabolomix.arf.tools import arf_plot_group_intensity
    session_state.session = session_state.AnalysisSession()
    out = json.loads(arf_plot_group_intensity(items=["PG"], groups=["ctrl"], output="payload"))
    assert out["error"]["code"] == "missing_state"


def test_payload_classes_names_and_nd(loaded):
    p = _payload(items=["PG", "PC", "PE", "PG 34:1"], groups=["ctrl", "ko"])
    by = {it["item"]: it for it in p["items"]}
    assert [s["spot_id"] for s in by["PG"]["spots"]] == [0, 1, 4]
    assert [s["spot_id"] for s in by["PC"]["spots"]] == [3]          # 内部標準は除く
    assert by["PC"]["n_spots_msms"] == 0
    assert by["PE"]["detected"] is False
    assert [s["spot_id"] for s in by["PG 34:1"]["spots"]] == [0]
    assert [g["label"] for g in p["groups"]] == ["ctrl", "ko"]
    assert p["groups"][0]["samples"] == ["20261006_ctrl_1", "20261006_ctrl_2"]


def test_standard_samples_drop_standard_only_spots(loaded):
    p = _payload(items=["PG"], groups=["ctrl", "ko"], standard_samples=["std"])
    assert [s["spot_id"] for s in p["items"][0]["spots"]] == [0, 1]
    assert p["excluded"]["standard_only"] == ["#4 PG 31:1|PG 17:0_14:1"]


def test_curation_flags_are_applied_by_default(loaded):
    store = curation_flags.FlagStore(curation_flags.curation_dir(loaded["arf2"]))
    key = curation_flags.alignment_key(loaded["arf2"])
    store.append([{"spot_id": 1, "flag": "wrong", "note": ""}], alignment=key,
                 review_id=None, source="test")
    p = _payload(items=["PG"], groups=["ctrl"])
    assert 1 not in [s["spot_id"] for s in p["items"][0]["spots"]]
    p = _payload(items=["PG"], groups=["ctrl"], apply_curation=False)
    assert 1 in [s["spot_id"] for s in p["items"][0]["spots"]]


def test_standard_filter_still_works_when_standard_sample_is_arf_excluded(loaded):
    session_state.session.arf.excluded_samples = {"20261006_std_1"}
    p = _payload(items=["PG"], groups=["ctrl", "ko"], standard_samples=["std"])
    assert [s["spot_id"] for s in p["items"][0]["spots"]] == [0, 1]
    assert p["excluded"]["standard_only"] == ["#4 PG 31:1|PG 17:0_14:1"]


def _write_review(arf2, review_id, sha, spots):
    directory = curation_flags.curation_dir(arf2)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"review-{review_id}.json").write_text(
        json.dumps({"alignment": {"alignment_sha256": sha}, "spots": spots}), encoding="utf-8")


def test_auto_likely_wrong_with_review_excludes_spots_and_ignores_other_alignment(loaded):
    arf2 = loaded["arf2"]
    sha = curation_flags.alignment_key(arf2)["alignment_sha256"]
    _write_review(arf2, "cr-20261006-120000-aaaa", sha, [{"spot_id": 1, "verdict": "likely_wrong"}])
    # 別アラインメント（sha256 違い）のレビューは、より新しくても無視される
    _write_review(arf2, "cr-20261006-130000-bbbb", "0" * 64, [{"spot_id": 0, "verdict": "likely_wrong"}])
    p = _payload(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert [s["spot_id"] for s in p["items"][0]["spots"]] == [0, 4]
    assert [e for e in p["excluded"]["auto_likely_wrong"]][0].startswith("#1 ")
    assert not any("レビュー" in c for c in p["caveats"])
    p = _payload(items=["PG"], groups=["ctrl"])
    assert 1 in [s["spot_id"] for s in p["items"][0]["spots"]]


def test_broken_flag_file_returns_error_with_message(loaded):
    directory = curation_flags.curation_dir(loaded["arf2"])
    directory.mkdir(parents=True, exist_ok=True)
    (directory / curation_flags.FLAGS_FILENAME).write_text("{not json\n", encoding="utf-8")
    out = _payload(items=["PG"], groups=["ctrl"])
    assert out["status"] == "error" and out["message"]


def test_auto_likely_wrong_without_review_adds_caveat(loaded):
    p = _payload(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert any("レビュー" in c for c in p["caveats"])


def test_excluded_samples_and_low_reliability(loaded):
    session_state.session.arf.excluded_samples = {"20261006_ko_2"}
    p = _payload(items=["PG"], groups=["ctrl", "ko"], low_reliability_samples=["ctrl_2"])
    assert p["groups"][1]["samples"] == ["20261006_ko_1"]
    ctrl = p["items"][0]["groups"][0]
    assert ctrl["n_in_stats"] == 1 and ctrl["samples"][1]["low_reliability"] is True


def test_group_fully_excluded_returns_error_json(loaded):
    session_state.session.arf.excluded_samples = {"20261006_ko_1", "20261006_ko_2"}
    out = _payload(items=["PG"], groups=["ctrl", "ko"])
    assert out["status"] == "error" and "No sample matched spec 'ko'" in out["message"]


def test_excluded_spots_are_not_drawn_in_class_or_name_parts(loaded):
    session_state.session.arf.excluded_spots = {1}
    p = _payload(items=["PG", "PG 35:1"], groups=["ctrl"])
    by = {it["item"]: it for it in p["items"]}
    assert 1 not in [s["spot_id"] for s in by["PG"]["spots"]]
    assert by["PG 35:1"]["spots"] == [] and by["PG 35:1"]["detected"] is False
    assert p["excluded"]["manual"] == ["#1 PG 35:1|PG 16:0_19:1"]


def test_caption_reports_manual_exclusions(loaded):
    from metabolomix.arf.tools import arf_plot_group_intensity
    session_state.session.arf.excluded_spots = {1}
    out = arf_plot_group_intensity(items=["PG"], groups=["ctrl"])
    assert "手動除外 1" in out[0]


def test_failed_call_discards_the_previous_figure(loaded):
    from metabolomix.arf.tools import arf_plot_group_intensity
    from metabolomix.tools.reports import save_figure
    arf_plot_group_intensity(items=["PG"], groups=["ctrl"], output="payload")
    assert session_state.session.arf.last_group_intensity is not None
    bad = _payload(items=["PG"], groups=["nonexistent"])
    assert bad["status"] == "error"
    assert session_state.session.arf.last_group_intensity is None
    out = json.loads(save_figure("group_intensity", "x"))
    assert out["error"]["code"] == "missing_state"


@pytest.mark.parametrize("kw", [{"detection_limit": -5}, {"detection_limit": 0},
                                {"detection_limit": float("nan")}, {"detection_limit": float("inf")},
                                {"ncols": 0}, {"ncols": 1.5}])
def test_invalid_detection_limit_or_ncols_returns_error_and_clears_slot(loaded, kw):
    from metabolomix.arf.tools import arf_plot_group_intensity
    arf_plot_group_intensity(items=["PG"], groups=["ctrl"], output="payload")
    out = _payload(items=["PG"], groups=["ctrl"], **kw)
    assert out["status"] == "error" and out["message"]
    assert session_state.session.arf.last_group_intensity is None


def test_image_mode_does_not_store_the_figure_when_render_fails(loaded, monkeypatch):
    from metabolomix.arf.tools import arf_plot_group_intensity
    from metabolomix.plots import group_intensity as gi

    def boom(*a, **k):
        raise RuntimeError("render failed")
    monkeypatch.setattr(gi, "render_group_intensity_plot", boom)
    with pytest.raises(RuntimeError):
        arf_plot_group_intensity(items=["PG"], groups=["ctrl"])
    assert session_state.session.arf.last_group_intensity is None


def test_invalid_output_returns_error_json(loaded):
    from metabolomix.arf.tools import arf_plot_group_intensity
    out = json.loads(arf_plot_group_intensity(items=["PG"], groups=["ctrl"], output="svg"))
    assert out["status"] == "error" and "output" in out["message"]


def test_detection_limit_from_param_file_and_argument(loaded, tmp_path):
    (tmp_path / "Dataset_x_param_202610061200.txt").write_text(
        "Minimum peak height: 1000\n", encoding="utf-8")
    p = _payload(items=["PG"], groups=["ctrl"])
    assert p["detection_limit"] == 1000.0 and p["detection_limit_source"] == "param_file"
    p = _payload(items=["PG"], groups=["ctrl"], detection_limit=500)
    assert p["detection_limit"] == 500.0 and p["detection_limit_source"] == "argument"


def test_image_output_and_session_keeps_payload(loaded):
    from mcp.server.fastmcp import Image
    from metabolomix.arf.tools import arf_plot_group_intensity
    out = arf_plot_group_intensity(items=["PG", "PE"], groups=["ctrl", "ko"])
    assert isinstance(out, list) and isinstance(out[1], Image)
    assert "2 項目" in out[0] and "N.D.: PE" in out[0]
    assert "PG 3 (MS/MS 0)" in out[0] and "PE 0" in out[0]
    assert session_state.session.arf.last_group_intensity["payload"]["items"][0]["item"] == "PG"


def test_save_without_plot_is_missing_state():
    from metabolomix.tools.reports import save_figure
    session_state.session = session_state.AnalysisSession()
    out = json.loads(save_figure("group_intensity", "x"))
    assert out["error"]["required_tools"] == ["arf_plot_group_intensity"]


def test_save_writes_png_and_svg(loaded, tmp_path, monkeypatch):
    from metabolomix.arf.tools import arf_plot_group_intensity
    from metabolomix.tools.reports import save_figure
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    arf_plot_group_intensity(items=["PG"], groups=["ctrl", "ko"], output="payload")
    msg = save_figure("group_intensity", "EV membrane")
    figures = tmp_path / "reports" / "figures"
    pngs = list(figures.glob("*_group_intensity.png"))
    assert len(pngs) == 1 and pngs[0].with_suffix(".svg").is_file()
    assert "group_intensity.png" in msg


def test_auto_likely_wrong_keeps_confirmed_spots(loaded):
    from metabolomix.msdial import tags as msdial_tags
    arf2 = loaded["arf2"]
    sha = curation_flags.alignment_key(arf2)["alignment_sha256"]
    _write_review(arf2, "cr-20261006-120000-aaaa", sha, [{"spot_id": 1, "verdict": "likely_wrong"}])
    msdial_tags.update_alignment_tag(msdial_tags.alignment_tag_path(arf2),
                                     tag_id=msdial_tags.CONFIRMED_TAG_ID, add=[1], remove=[])
    p = _payload(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert 1 in [s["spot_id"] for s in p["items"][0]["spots"]]
    assert p["excluded"]["auto_likely_wrong"] == []
