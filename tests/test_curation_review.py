import json
import shutil
import subprocess

import pytest

from lipidmix.arf2.reader import load_catalog
from lipidmix.curation import evidence, flags, judge, review, viewer
from lipidmix.library import store as library_store
from tests.curation_fixtures import write_alignment_set


@pytest.fixture()
def built(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg")
    s = library_store.open_store(paths["msp"])
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    result = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                               th=judge.resolve_thresholds(None), file_ids=None, max_traces=12,
                               selection={"kind": "annotated"})
    yield paths, result
    s.close()


def test_verdicts_follow_the_evidence(built):
    _, result = built
    by_id = {s["spot_id"]: s for s in result["spots"]}
    assert by_id[0]["verdict"] == "ok"
    # spot 1: low_score（弱）+ drt_out 1.5 分（弱）+ ギザギザの EIC（弱）→ suspect
    assert by_id[1]["verdict"] == "suspect"
    assert {"low_score", "drt_out", "eic_poor"} <= set(by_id[1]["reasons"])
    assert result["counts"] == {"ok": 1, "suspect": 1, "likely_wrong": 0}


def test_summary_tsv_lists_only_flagged_spots(built):
    _, result = built
    lines = review.summary_tsv(result).splitlines()
    assert lines[0].split("\t")[:3] == ["spot_id", "name", "ontology"]
    assert [line.split("\t")[0] for line in lines[1:]] == ["1"]


def test_save_and_load_round_trip(built):
    paths, result = built
    saved = review.save_review(result)
    assert saved["json"].parent == flags.curation_dir(paths["arf2"])
    loaded = review.load_review(paths["arf2"], result["review_id"])
    assert loaded["review_id"] == result["review_id"]
    html = saved["html"].read_text(encoding="utf-8")
    assert "</script>" in html
    assert result["review_id"] in html


def test_embedded_json_cannot_break_out_of_the_script_tag(built):
    _, result = built
    result["spots"][0]["name"] = "</script><b>x</b>"
    html = viewer.render_html(result)
    assert "</script><b>" not in html


def test_existing_flags_are_attached_to_spots(built):
    paths, result = built
    store = flags.FlagStore(flags.curation_dir(paths["arf2"]))
    store.append([{"spot_id": 1, "flag": "wrong"}], alignment=result["alignment"],
                 review_id="old", source="user")
    s = library_store.open_store(paths["msp"])
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        again = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                  th=judge.resolve_thresholds(None), file_ids=None,
                                  max_traces=12, selection={"kind": "annotated"})
    finally:
        s.close()
    assert {sp["spot_id"]: sp["flag"] for sp in again["spots"]} == {0: None, 1: "wrong"}


def test_low_reference_resolution_warns_about_the_library(built, tmp_path):
    paths, _ = built
    other = tmp_path / "other.msp"
    other.write_text("NAME: X\nPRECURSORMZ: 100.0\nIONMODE: Positive\nNum Peaks: 0\n", encoding="utf-8")
    s = library_store.open_store(other)
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        result = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                   th=judge.resolve_thresholds(None), file_ids=None,
                                   max_traces=12, selection={"kind": "annotated"})
    finally:
        s.close()
    assert any("別のライブラリ" in w for w in result["warnings"])


def test_page_splits_spots(built):
    _, result = built
    first = review.page(result, 0, page_size=1)
    assert first["n_pages"] == 2 and len(first["spots"]) == 1


def test_app_template_has_no_embedded_data():
    html = viewer.render_html(None)
    assert "const EMBEDDED = null" in html


def test_template_avoids_horizontal_overflow_at_phone_width():
    # 実ファイル名にはスペースがなく、header h1 が縮められないと sticky header ごと
    # 横スクロールが出る（binding requirement: no horizontal page scroll at phone
    # width）。#trends / #grid の minmax も 320px 幅の携帯（16px ガター）では
    # 300/320px 未満に縮まないと同様に溢れる。
    html = viewer.render_html(None)
    assert "overflow-wrap:anywhere" in html
    assert "min(300px, 100%)" in html


def test_render_builds_cards_in_chunks_with_a_cancellable_generation_counter():
    # 1468 件のような大きなレビューで render() が 1 回の同期ループで全カードを
    # 組み立てると、ページが 10-15 秒応答しなくなる（実データ check）。チャンク分割し、
    # フィルタ変更などで新しい render() が割り込んだら古い方を打ち切る世代カウンタが要る。
    html = viewer.render_html(None)
    assert "renderGeneration" in html
    assert "generation !== renderGeneration" in html
    assert "requestAnimationFrame" in html
    assert "setTimeout" in html
    assert "DocumentFragment" in html


def test_draw_trend_clips_the_regression_line_to_the_plot_rectangle():
    # slope*x+intercept は RT レンジの外まで伸びうり、パネルからはみ出して見える
    # （実データ check）。プロット矩形で ctx.clip() する。
    html = viewer.render_html(None)
    start = html.index("function drawTrend")
    end = html.index("function spotCard")
    trend_source = html[start:end]
    assert "ctx.save()" in trend_source
    assert "ctx.clip()" in trend_source
    assert "ctx.restore()" in trend_source


def test_extracted_script_is_valid_javascript(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が見つからないので構文チェックを省略")
    html = viewer.render_html(None)
    start = html.index("<script>") + len("<script>")
    end = html.index("</script>")
    script_path = tmp_path / "viewer.js"
    script_path.write_text(html[start:end], encoding="utf-8")
    result = subprocess.run([node, "--check", str(script_path)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# ---------- I2: TSV の行数上限 ----------

def test_summary_tsv_keeps_the_first_max_rows_by_severity(built):
    _, result = built
    result["spots"][0]["flag"] = "suspect"      # ok だがフラグ済み → 表に載る(2 行になる)
    assert review.n_summary_rows(result) == 2
    lines = review.summary_tsv(result, max_rows=1).splitlines()
    assert len(lines) == 2                      # 見出し + 1 行
    assert lines[1].split("\t")[0] == "1"       # suspect が先(重い順)


# ---------- M1: review_id の形式 ----------

@pytest.mark.parametrize("bad", ["../x", "cr-1", "cr-20260101-000000-zzzz", "cr-20260101-000000-abcd/..", ""])
def test_load_review_rejects_malformed_review_ids(built, bad):
    paths, _ = built
    with pytest.raises(ValueError):
        review.load_review(paths["arf2"], bad)


def test_generated_review_ids_are_valid(built):
    _, result = built
    assert review.is_valid_review_id(result["review_id"])


# ---------- M5: 置換の順序と `<` のエスケープ ----------

def test_spot_names_containing_the_prefix_placeholder_survive(built):
    _, result = built
    result["spots"][0]["name"] = "__SUBMISSION_PREFIX__ <!-- x"
    html = viewer.render_html(result)
    start = html.index("const EMBEDDED = ") + len("const EMBEDDED = ")
    end = html.index(";\nconst PREFIX")
    embedded = html[start:end]
    assert "<" not in embedded
    assert json.loads(embedded)["spots"][0]["name"] == "__SUBMISSION_PREFIX__ <!-- x"


# ---------- M7: クラス別傾向の要約 ----------

def test_trend_summary_includes_outlier_counts():
    fake = {"trend": {"classes": {"PC": {"n": 9, "r2": 0.95, "n_outliers": 1, "coef": [1, 2]}}}}
    assert review.trend_summary(fake) == {"PC": {"n": 9, "r2": 0.95, "n_outliers": 1}}


# ---------- I4 / I5 / M4: ビューアのテンプレート ----------

def _script():
    html = viewer.render_html(None)
    return html[html.index("<script>"):html.index("</script>")]


def test_submission_includes_the_arf2_path():
    script = _script()
    body = script[script.index("function submission"):]
    body = body[:body.index("\n}")]
    assert "REVIEW.arf2_path" in body


def test_trend_canvas_click_jumps_to_the_spot_card():
    script = _script()
    assert "jumpToSpot" in script
    trends = script[script.index("function renderTrends"):]
    trends = trends[:trends.index("\n}")]
    assert 'addEventListener("click"' in trends
    assert "HIT_RADIUS" in script
    jump = script[script.index("function jumpToSpot"):]
    jump = jump[:jump.index("\n}")]
    assert "scrollIntoView" in jump and "highlight" in jump
    assert "resetFilters" in jump


def test_card_canvases_open_an_enlarged_dialog():
    html = viewer.render_html(None)
    assert '<dialog id="zoom"' in html
    script = _script()
    card = script[script.index("function spotCard"):]
    card = card[:card.index("\n}")]
    assert "openZoom" in card
    zoom = script[script.index("function openZoom"):]
    zoom = zoom[:zoom.index("\n}")]
    assert "drawEic" in zoom and "drawMirror" in zoom
    assert "Escape" in script
    assert 'id="zoom-close"' in html


def test_note_input_reflects_unsent_edits():
    script = _script()
    card = script[script.index("function spotCard"):]
    card = card[:card.index("\n}")]
    assert "edits.get(spot.spot_id).note" in card
    assert "spot.mz?.toFixed" not in card


def test_app_loader_shows_error_payloads_and_send_updates_spots():
    script = _script()
    assert "showError" in script
    send = script[script.index('getElementById("send").addEventListener'):]
    send = send[:send.index("\n});")]
    assert "spot.flag = " in send and "spot.flag_note = " in send
    assert send.index("spot.flag = ") < send.index("edits.clear()")
